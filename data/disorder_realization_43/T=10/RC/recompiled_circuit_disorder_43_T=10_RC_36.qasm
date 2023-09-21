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
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1057518) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(-2.7266154) q[0];
x q[1];
rz(0.33306723) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.2531467) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9589899) q[1];
sx q[1];
rz(-0.62287736) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(1.9777771) q[3];
sx q[3];
rz(-2.3325936) q[3];
sx q[3];
rz(0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(0.17856199) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-0.006342412) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(2.1727174) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(2.7768163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9176863) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(2.9427337) q[1];
rz(-pi) q[2];
rz(0.56224058) q[3];
sx q[3];
rz(-1.3573109) q[3];
sx q[3];
rz(-1.7052887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851615) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(-3.1306144) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(0.095741622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75671065) q[0];
sx q[0];
rz(-1.7236992) q[0];
sx q[0];
rz(-1.3923313) q[0];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.1635457) q[2];
sx q[2];
rz(0.72088748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0248191) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(-2.5953672) q[1];
rz(-pi) q[2];
rz(0.057283244) q[3];
sx q[3];
rz(-0.55192845) q[3];
sx q[3];
rz(-0.79997593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(2.1077164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(0.69181504) q[0];
x q[1];
rz(-2.6150377) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(-1.9499792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.911072) q[1];
sx q[1];
rz(-1.5057179) q[1];
sx q[1];
rz(1.7622403) q[1];
x q[2];
rz(1.2437808) q[3];
sx q[3];
rz(-1.838284) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(2.4436061) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-0.57410747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363279) q[0];
sx q[0];
rz(-2.5308373) q[0];
sx q[0];
rz(1.3461793) q[0];
rz(-pi) q[1];
rz(-1.5482043) q[2];
sx q[2];
rz(-0.85561692) q[2];
sx q[2];
rz(-1.9695645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9088604) q[1];
sx q[1];
rz(-0.4545916) q[1];
sx q[1];
rz(-2.4652387) q[1];
x q[2];
rz(-1.8168713) q[3];
sx q[3];
rz(-2.4613791) q[3];
sx q[3];
rz(2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(2.4093157) q[1];
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
x q[1];
rz(0.90768355) q[2];
sx q[2];
rz(-2.216202) q[2];
sx q[2];
rz(-1.6931319) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9975035) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-1.3527649) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(0.41527173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419287) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(1.3616189) q[0];
rz(-pi) q[1];
rz(-0.98203512) q[2];
sx q[2];
rz(-1.2274449) q[2];
sx q[2];
rz(-2.8729168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60939497) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(2.1983912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9023864) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(2.1945206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-2.1264123) q[2];
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
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-0.24169895) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1358571) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(0.70941305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62944062) q[2];
sx q[2];
rz(-1.9295613) q[2];
sx q[2];
rz(1.7164873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.3543345) q[1];
x q[2];
rz(1.6493158) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1239132) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(0.29385847) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732849) q[0];
sx q[0];
rz(-2.377016) q[0];
sx q[0];
rz(1.1029878) q[0];
x q[1];
rz(0.80957885) q[2];
sx q[2];
rz(-2.5123345) q[2];
sx q[2];
rz(1.6639683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2453354) q[1];
sx q[1];
rz(-0.63071139) q[1];
sx q[1];
rz(-2.3920822) q[1];
rz(-pi) q[2];
rz(-2.1095554) q[3];
sx q[3];
rz(-1.7146535) q[3];
sx q[3];
rz(1.1080351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(0.59797257) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(1.3806608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1191694) q[0];
sx q[0];
rz(-1.701231) q[0];
sx q[0];
rz(-2.2020257) q[0];
rz(-pi) q[1];
rz(-1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.6053111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8163562) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(-0.34646323) q[1];
rz(-pi) q[2];
rz(1.7396183) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(-1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(-2.5972988) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(-2.9059698) q[3];
sx q[3];
rz(-1.0141254) q[3];
sx q[3];
rz(0.49401382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
