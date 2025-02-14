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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(-0.21633202) q[0];
rz(4.0989838) q[1];
sx q[1];
rz(5.6070072) q[1];
sx q[1];
rz(11.054872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0441598) q[0];
sx q[0];
rz(-2.5852647) q[0];
sx q[0];
rz(-0.018228368) q[0];
x q[1];
rz(1.2443107) q[2];
sx q[2];
rz(-2.0225581) q[2];
sx q[2];
rz(-1.8725852) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2547467) q[1];
sx q[1];
rz(-2.0664082) q[1];
sx q[1];
rz(1.1911456) q[1];
x q[2];
rz(2.7483398) q[3];
sx q[3];
rz(-2.234708) q[3];
sx q[3];
rz(2.8877986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1067074) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(-2.9599221) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492391) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-0.44678584) q[0];
rz(-2.2611639) q[1];
sx q[1];
rz(-1.3648405) q[1];
sx q[1];
rz(0.78278881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7004958) q[0];
sx q[0];
rz(-1.6320328) q[0];
sx q[0];
rz(1.815531) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5704284) q[2];
sx q[2];
rz(-0.70020247) q[2];
sx q[2];
rz(-0.68298662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0228717) q[1];
sx q[1];
rz(-2.0346271) q[1];
sx q[1];
rz(1.9074341) q[1];
x q[2];
rz(0.47329013) q[3];
sx q[3];
rz(-0.9447228) q[3];
sx q[3];
rz(-1.9260709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.030152628) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(-2.1622369) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6556743) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(-2.3517877) q[0];
rz(2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(-2.1479215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8756722) q[0];
sx q[0];
rz(-1.3215995) q[0];
sx q[0];
rz(2.5795291) q[0];
rz(1.7970071) q[2];
sx q[2];
rz(-0.96983428) q[2];
sx q[2];
rz(0.41659875) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34383306) q[1];
sx q[1];
rz(-2.0270837) q[1];
sx q[1];
rz(-0.47618687) q[1];
rz(-pi) q[2];
rz(2.42071) q[3];
sx q[3];
rz(-1.3647121) q[3];
sx q[3];
rz(-0.64180798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32187244) q[2];
sx q[2];
rz(-2.0237782) q[2];
sx q[2];
rz(0.37128386) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.4042847) q[0];
rz(0.67717254) q[1];
sx q[1];
rz(-1.9858457) q[1];
sx q[1];
rz(-1.4926532) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9335404) q[0];
sx q[0];
rz(-1.786199) q[0];
sx q[0];
rz(2.893704) q[0];
x q[1];
rz(2.8437496) q[2];
sx q[2];
rz(-1.9854416) q[2];
sx q[2];
rz(0.21098247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74977743) q[1];
sx q[1];
rz(-1.3484553) q[1];
sx q[1];
rz(-0.27395497) q[1];
x q[2];
rz(-1.221232) q[3];
sx q[3];
rz(-2.3858983) q[3];
sx q[3];
rz(0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(0.34823927) q[2];
rz(-1.6866775) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763879) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(-0.46420321) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(1.8468599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(-2.110092) q[0];
rz(-pi) q[1];
rz(1.4424397) q[2];
sx q[2];
rz(-2.1048628) q[2];
sx q[2];
rz(2.4160224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.059493493) q[1];
sx q[1];
rz(-1.1319185) q[1];
sx q[1];
rz(-2.782269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5312349) q[3];
sx q[3];
rz(-0.98516432) q[3];
sx q[3];
rz(-2.741339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(-0.56582212) q[2];
rz(-0.081341751) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(0.62059039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6767947) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(1.5860522) q[0];
rz(1.0606891) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(2.5097844) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.466296) q[0];
sx q[0];
rz(-1.6501556) q[0];
sx q[0];
rz(2.9064889) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.00497) q[2];
sx q[2];
rz(-2.3477051) q[2];
sx q[2];
rz(-2.6086406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28367701) q[1];
sx q[1];
rz(-1.9092968) q[1];
sx q[1];
rz(-1.0033016) q[1];
rz(-2.9618044) q[3];
sx q[3];
rz(-2.4500812) q[3];
sx q[3];
rz(2.9464242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98171988) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(1.8286797) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(-1.7074728) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(-0.69951192) q[0];
rz(-0.46547678) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(-2.2043998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6420744) q[0];
sx q[0];
rz(-2.8694911) q[0];
sx q[0];
rz(2.4689552) q[0];
rz(-pi) q[1];
rz(0.40528751) q[2];
sx q[2];
rz(-2.8817085) q[2];
sx q[2];
rz(1.056162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7286789) q[1];
sx q[1];
rz(-1.2866486) q[1];
sx q[1];
rz(0.22032622) q[1];
x q[2];
rz(-2.8392467) q[3];
sx q[3];
rz(-1.28966) q[3];
sx q[3];
rz(0.1880364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-0.74791351) q[0];
rz(-0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(-0.0029729923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2073313) q[0];
sx q[0];
rz(-2.2483279) q[0];
sx q[0];
rz(0.24768655) q[0];
rz(-pi) q[1];
rz(-0.019549088) q[2];
sx q[2];
rz(-1.6416993) q[2];
sx q[2];
rz(-2.8343458) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5103448) q[1];
sx q[1];
rz(-1.6142134) q[1];
sx q[1];
rz(-1.7249291) q[1];
rz(-pi) q[2];
rz(1.3722754) q[3];
sx q[3];
rz(-1.7619507) q[3];
sx q[3];
rz(1.6894345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(-0.083960697) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(2.3626204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(0.10636605) q[0];
rz(0.97995177) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(0.76593691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1003215) q[0];
sx q[0];
rz(-2.092397) q[0];
sx q[0];
rz(-2.0727511) q[0];
rz(0.19975234) q[2];
sx q[2];
rz(-1.6589266) q[2];
sx q[2];
rz(2.5823808) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3043266) q[1];
sx q[1];
rz(-0.89511739) q[1];
sx q[1];
rz(1.6153687) q[1];
rz(2.8838653) q[3];
sx q[3];
rz(-3.1270087) q[3];
sx q[3];
rz(-0.38250438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6672259) q[2];
sx q[2];
rz(-2.6211278) q[2];
sx q[2];
rz(-2.4412947) q[2];
rz(-0.66655603) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(2.0535645) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67548442) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(1.745537) q[0];
rz(-1.4736157) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(-2.4748763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5294801) q[0];
sx q[0];
rz(-0.94769883) q[0];
sx q[0];
rz(-2.3701131) q[0];
x q[1];
rz(-0.36138968) q[2];
sx q[2];
rz(-2.4494684) q[2];
sx q[2];
rz(-0.036594242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4426431) q[1];
sx q[1];
rz(-1.2546228) q[1];
sx q[1];
rz(1.433028) q[1];
rz(-pi) q[2];
rz(2.1176841) q[3];
sx q[3];
rz(-1.0677665) q[3];
sx q[3];
rz(-1.5707113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0746158) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(-0.43205076) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(-2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074832) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(-1.8941849) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(0.16130372) q[2];
sx q[2];
rz(-1.9111173) q[2];
sx q[2];
rz(2.5572122) q[2];
rz(1.5374938) q[3];
sx q[3];
rz(-1.0819482) q[3];
sx q[3];
rz(-0.29500189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
