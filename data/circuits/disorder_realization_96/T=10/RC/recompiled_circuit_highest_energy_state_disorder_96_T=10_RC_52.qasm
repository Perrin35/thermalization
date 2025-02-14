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
rz(0.32956707) q[0];
sx q[0];
rz(-2.1962533) q[0];
sx q[0];
rz(-3.1402631) q[0];
rz(0.20345649) q[1];
sx q[1];
rz(-0.23509547) q[1];
sx q[1];
rz(-2.5703365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205367) q[0];
sx q[0];
rz(-2.1696496) q[0];
sx q[0];
rz(-2.1355892) q[0];
x q[1];
rz(-0.18988256) q[2];
sx q[2];
rz(-1.7203091) q[2];
sx q[2];
rz(0.25611311) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3498342) q[1];
sx q[1];
rz(-1.5538284) q[1];
sx q[1];
rz(-2.6512673) q[1];
x q[2];
rz(1.9523078) q[3];
sx q[3];
rz(-1.1433868) q[3];
sx q[3];
rz(1.007146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0433537) q[2];
sx q[2];
rz(-1.9184155) q[2];
sx q[2];
rz(-2.5377048) q[2];
rz(-1.0107001) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(-2.219626) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-1.0930467) q[0];
sx q[0];
rz(0.0086722886) q[0];
rz(1.1588833) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(-0.11223665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771292) q[0];
sx q[0];
rz(-2.5647129) q[0];
sx q[0];
rz(0.5441027) q[0];
rz(-pi) q[1];
rz(2.49493) q[2];
sx q[2];
rz(-1.5407667) q[2];
sx q[2];
rz(0.5867556) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48277897) q[1];
sx q[1];
rz(-2.5943752) q[1];
sx q[1];
rz(-2.7040255) q[1];
rz(-2.3149764) q[3];
sx q[3];
rz(-0.38020779) q[3];
sx q[3];
rz(2.2728863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9021641) q[2];
sx q[2];
rz(-1.8205234) q[2];
sx q[2];
rz(0.93977896) q[2];
rz(2.2043665) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(-2.7711813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060191) q[0];
sx q[0];
rz(-2.2414099) q[0];
sx q[0];
rz(-2.5584333) q[0];
rz(1.0519625) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(-3.1254056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.095854) q[0];
sx q[0];
rz(-1.7940137) q[0];
sx q[0];
rz(2.2824817) q[0];
rz(-2.8444958) q[2];
sx q[2];
rz(-1.8645745) q[2];
sx q[2];
rz(-2.370594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0250108) q[1];
sx q[1];
rz(-1.407928) q[1];
sx q[1];
rz(3.0322269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82456581) q[3];
sx q[3];
rz(-2.0494132) q[3];
sx q[3];
rz(-2.5156817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6215324) q[2];
sx q[2];
rz(-1.748087) q[2];
sx q[2];
rz(0.072546093) q[2];
rz(-1.0091311) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(-0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9811454) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(-2.2571046) q[0];
rz(1.6403987) q[1];
sx q[1];
rz(-1.4720935) q[1];
sx q[1];
rz(-2.5515058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9035066) q[0];
sx q[0];
rz(-1.0732539) q[0];
sx q[0];
rz(2.7142798) q[0];
x q[1];
rz(-0.22802148) q[2];
sx q[2];
rz(-1.8393901) q[2];
sx q[2];
rz(3.0326469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9635439) q[1];
sx q[1];
rz(-2.5317988) q[1];
sx q[1];
rz(-2.1765373) q[1];
x q[2];
rz(-2.8899412) q[3];
sx q[3];
rz(-2.6918937) q[3];
sx q[3];
rz(0.77317274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13714743) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(-1.2443292) q[2];
rz(1.7826049) q[3];
sx q[3];
rz(-0.46670306) q[3];
sx q[3];
rz(0.14314237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2537848) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(-2.7746871) q[0];
rz(1.3430355) q[1];
sx q[1];
rz(-0.29452205) q[1];
sx q[1];
rz(1.8273182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2283233) q[0];
sx q[0];
rz(-2.0073038) q[0];
sx q[0];
rz(-0.69605791) q[0];
x q[1];
rz(1.8601244) q[2];
sx q[2];
rz(-2.0380033) q[2];
sx q[2];
rz(1.510646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6180743) q[1];
sx q[1];
rz(-0.95361751) q[1];
sx q[1];
rz(0.38477044) q[1];
x q[2];
rz(-1.9461104) q[3];
sx q[3];
rz(-1.2367782) q[3];
sx q[3];
rz(0.53900146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5985976) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(-2.561595) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5744659) q[3];
sx q[3];
rz(-2.4549761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065303236) q[0];
sx q[0];
rz(-2.6798601) q[0];
sx q[0];
rz(2.8522458) q[0];
rz(-3.0416327) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(0.79773703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432309) q[0];
sx q[0];
rz(-1.7055178) q[0];
sx q[0];
rz(-0.93905296) q[0];
rz(-2.9399273) q[2];
sx q[2];
rz(-1.4603134) q[2];
sx q[2];
rz(-0.1486272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9790837) q[1];
sx q[1];
rz(-1.6220656) q[1];
sx q[1];
rz(-2.8460345) q[1];
rz(-pi) q[2];
rz(-2.3008505) q[3];
sx q[3];
rz(-1.4952097) q[3];
sx q[3];
rz(-0.10160916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37819698) q[2];
sx q[2];
rz(-1.8969456) q[2];
sx q[2];
rz(2.3217616) q[2];
rz(1.648692) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(2.7259887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019526871) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(2.8269738) q[0];
rz(0.64612359) q[1];
sx q[1];
rz(-1.6982634) q[1];
sx q[1];
rz(0.12166469) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3356161) q[0];
sx q[0];
rz(-1.8467963) q[0];
sx q[0];
rz(-1.9673303) q[0];
rz(-pi) q[1];
rz(-0.88603772) q[2];
sx q[2];
rz(-2.2412063) q[2];
sx q[2];
rz(-1.3051906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77326425) q[1];
sx q[1];
rz(-0.46302893) q[1];
sx q[1];
rz(-2.8431588) q[1];
rz(-pi) q[2];
rz(0.1945576) q[3];
sx q[3];
rz(-2.1702507) q[3];
sx q[3];
rz(0.50473467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.17522) q[2];
sx q[2];
rz(-0.97869527) q[2];
rz(0.77224246) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(-2.1286807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.1145645) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(1.8373328) q[0];
rz(-0.14106855) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(1.28349) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9314037) q[0];
sx q[0];
rz(-2.922466) q[0];
sx q[0];
rz(1.7396084) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2705994) q[2];
sx q[2];
rz(-1.6643447) q[2];
sx q[2];
rz(1.5644422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56399689) q[1];
sx q[1];
rz(-1.3118) q[1];
sx q[1];
rz(1.6902655) q[1];
x q[2];
rz(0.047824511) q[3];
sx q[3];
rz(-2.3615814) q[3];
sx q[3];
rz(-0.081950233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86113247) q[2];
sx q[2];
rz(-2.2304163) q[2];
sx q[2];
rz(-0.99315161) q[2];
rz(-1.447621) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(-1.2710458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3250378) q[0];
sx q[0];
rz(-3.0386381) q[0];
sx q[0];
rz(-2.2651941) q[0];
rz(1.5429629) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(-2.0230944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3680822) q[0];
sx q[0];
rz(-2.4187669) q[0];
sx q[0];
rz(2.5834492) q[0];
rz(-pi) q[1];
rz(-0.72965168) q[2];
sx q[2];
rz(-2.0579946) q[2];
sx q[2];
rz(-1.2318357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1269605) q[1];
sx q[1];
rz(-1.4585969) q[1];
sx q[1];
rz(0.87027624) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7574945) q[3];
sx q[3];
rz(-2.4599061) q[3];
sx q[3];
rz(-3.0317999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2529651) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(1.2523119) q[2];
rz(-2.0400932) q[3];
sx q[3];
rz(-1.3809729) q[3];
sx q[3];
rz(1.8496877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0970704) q[0];
sx q[0];
rz(-3.0945393) q[0];
sx q[0];
rz(0.72730056) q[0];
rz(-0.76668382) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(-1.1904221) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.577548) q[0];
sx q[0];
rz(-2.2909031) q[0];
sx q[0];
rz(1.4029361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4486362) q[2];
sx q[2];
rz(-2.4776931) q[2];
sx q[2];
rz(0.66421504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1228635) q[1];
sx q[1];
rz(-1.3236538) q[1];
sx q[1];
rz(-0.31063947) q[1];
x q[2];
rz(-0.50325692) q[3];
sx q[3];
rz(-1.0531593) q[3];
sx q[3];
rz(2.518017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7850354) q[2];
sx q[2];
rz(-2.0544453) q[2];
sx q[2];
rz(-1.1615151) q[2];
rz(-2.0098497) q[3];
sx q[3];
rz(-2.3029906) q[3];
sx q[3];
rz(1.0034126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3840735) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(-1.1643428) q[2];
sx q[2];
rz(-0.60383527) q[2];
sx q[2];
rz(2.8964103) q[2];
rz(-1.499621) q[3];
sx q[3];
rz(-0.19656062) q[3];
sx q[3];
rz(1.82609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
