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
rz(1.0951618) q[0];
sx q[0];
rz(-2.996063) q[0];
sx q[0];
rz(-0.50330436) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(1.4454747) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88004011) q[0];
sx q[0];
rz(-1.6140811) q[0];
sx q[0];
rz(-0.39899428) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93930556) q[2];
sx q[2];
rz(-0.92512265) q[2];
sx q[2];
rz(-1.948057) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9564285) q[1];
sx q[1];
rz(-1.899029) q[1];
sx q[1];
rz(-2.4976455) q[1];
rz(-pi) q[2];
rz(1.5551994) q[3];
sx q[3];
rz(-1.6199779) q[3];
sx q[3];
rz(1.552207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4989) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(-1.5042245) q[2];
rz(0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40642834) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(-0.51097393) q[0];
rz(-1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-2.8151292) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24169479) q[0];
sx q[0];
rz(-2.7621671) q[0];
sx q[0];
rz(2.3510758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9316177) q[2];
sx q[2];
rz(-2.2433503) q[2];
sx q[2];
rz(-0.58174101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4703456) q[1];
sx q[1];
rz(-1.2675646) q[1];
sx q[1];
rz(0.95945759) q[1];
x q[2];
rz(2.191698) q[3];
sx q[3];
rz(-0.46484676) q[3];
sx q[3];
rz(-2.6863033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9713126) q[2];
sx q[2];
rz(-0.22394094) q[2];
sx q[2];
rz(-1.2962606) q[2];
rz(0.41804677) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068785) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(0.54247722) q[0];
rz(-1.7659278) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-3.1105522) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2707996) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(-0.55678456) q[0];
x q[1];
rz(2.228187) q[2];
sx q[2];
rz(-1.3904461) q[2];
sx q[2];
rz(3.0424398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2636288) q[1];
sx q[1];
rz(-2.1641556) q[1];
sx q[1];
rz(0.29747648) q[1];
rz(-0.026588566) q[3];
sx q[3];
rz(-0.55517653) q[3];
sx q[3];
rz(0.41586499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5151908) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(-2.314563) q[2];
rz(2.6229897) q[3];
sx q[3];
rz(-2.0293472) q[3];
sx q[3];
rz(1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14722918) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.2116785) q[0];
rz(-1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.3406219) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3163534) q[0];
sx q[0];
rz(-1.2636501) q[0];
sx q[0];
rz(0.65748837) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.421809) q[2];
sx q[2];
rz(-1.8049311) q[2];
sx q[2];
rz(-2.9643167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2681899) q[1];
sx q[1];
rz(-2.4854069) q[1];
sx q[1];
rz(2.6471958) q[1];
x q[2];
rz(-0.31130917) q[3];
sx q[3];
rz(-1.9640067) q[3];
sx q[3];
rz(-1.1953199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2010605) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(2.000467) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-0.96735668) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56213266) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(-0.60254565) q[0];
rz(2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(0.7811195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8666894) q[0];
sx q[0];
rz(-2.8056228) q[0];
sx q[0];
rz(2.1076403) q[0];
x q[1];
rz(0.94046142) q[2];
sx q[2];
rz(-2.2302719) q[2];
sx q[2];
rz(2.3465921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32051099) q[1];
sx q[1];
rz(-0.85144224) q[1];
sx q[1];
rz(-0.11493747) q[1];
rz(2.9424465) q[3];
sx q[3];
rz(-0.73561397) q[3];
sx q[3];
rz(0.74935952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1802804) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(2.9126634) q[2];
rz(1.3198352) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(0.84407097) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22885403) q[0];
sx q[0];
rz(-1.5541394) q[0];
sx q[0];
rz(1.1786906) q[0];
rz(-1.9121869) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(1.2961402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4018009) q[0];
sx q[0];
rz(-1.8515046) q[0];
sx q[0];
rz(0.14174353) q[0];
rz(-2.3713114) q[2];
sx q[2];
rz(-2.2438223) q[2];
sx q[2];
rz(-1.1550624) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78627045) q[1];
sx q[1];
rz(-0.67364022) q[1];
sx q[1];
rz(0.73553971) q[1];
rz(0.48719897) q[3];
sx q[3];
rz(-1.4811084) q[3];
sx q[3];
rz(-1.7023757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65427762) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(2.8793092) q[2];
rz(-2.8413963) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(-0.20849553) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7009785) q[0];
sx q[0];
rz(-2.746026) q[0];
sx q[0];
rz(1.3025008) q[0];
rz(2.473096) q[1];
sx q[1];
rz(-1.786247) q[1];
sx q[1];
rz(1.0350593) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9205017) q[0];
sx q[0];
rz(-1.1701487) q[0];
sx q[0];
rz(2.9095838) q[0];
rz(-pi) q[1];
rz(1.9403753) q[2];
sx q[2];
rz(-1.7453004) q[2];
sx q[2];
rz(0.39350393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3429195) q[1];
sx q[1];
rz(-1.2274776) q[1];
sx q[1];
rz(-3.0613106) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5104654) q[3];
sx q[3];
rz(-2.2911117) q[3];
sx q[3];
rz(-2.9583954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11789007) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(-1.4855851) q[2];
rz(0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(-2.7070467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3561803) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(-2.8670512) q[0];
rz(-1.9369594) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(2.5801632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954079) q[0];
sx q[0];
rz(-2.2642379) q[0];
sx q[0];
rz(-2.2827143) q[0];
rz(-pi) q[1];
rz(-1.2677293) q[2];
sx q[2];
rz(-2.2243739) q[2];
sx q[2];
rz(-2.4619964) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16374396) q[1];
sx q[1];
rz(-1.1625966) q[1];
sx q[1];
rz(-0.331649) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9188324) q[3];
sx q[3];
rz(-1.819918) q[3];
sx q[3];
rz(0.083569738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56208912) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(-1.7378463) q[2];
rz(-1.7956519) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(-2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76029921) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(-2.2316933) q[0];
rz(1.9445885) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(-2.2534456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4662381) q[0];
sx q[0];
rz(-0.14396891) q[0];
sx q[0];
rz(-0.86021249) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65983628) q[2];
sx q[2];
rz(-1.5232695) q[2];
sx q[2];
rz(2.8328676) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8348114) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(2.4624834) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2751332) q[3];
sx q[3];
rz(-1.5348892) q[3];
sx q[3];
rz(1.1334238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1278648) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(-2.0590651) q[2];
rz(0.23076375) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312663) q[0];
sx q[0];
rz(-2.3895097) q[0];
sx q[0];
rz(2.5600774) q[0];
rz(0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(1.2967348) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49823495) q[0];
sx q[0];
rz(-1.5603934) q[0];
sx q[0];
rz(-1.5525981) q[0];
rz(-pi) q[1];
rz(-0.34411547) q[2];
sx q[2];
rz(-1.5617579) q[2];
sx q[2];
rz(1.8923541) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.30137128) q[1];
sx q[1];
rz(-2.6379963) q[1];
sx q[1];
rz(-2.1238158) q[1];
rz(-0.48892816) q[3];
sx q[3];
rz(-2.4081514) q[3];
sx q[3];
rz(1.6176002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7221308) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(2.6922743) q[2];
rz(0.32014534) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(-1.8951353) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37353361) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(-1.7569348) q[1];
sx q[1];
rz(-0.23778267) q[1];
sx q[1];
rz(2.3429088) q[1];
rz(0.25945406) q[2];
sx q[2];
rz(-2.1726407) q[2];
sx q[2];
rz(3.0214027) q[2];
rz(0.89775916) q[3];
sx q[3];
rz(-2.5644292) q[3];
sx q[3];
rz(0.226365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
