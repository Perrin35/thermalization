OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(-1.5919332) q[1];
sx q[1];
rz(-3.0631493) q[1];
sx q[1];
rz(0.6426386) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064518236) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(-1.5383188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.379001) q[2];
sx q[2];
rz(-1.0704874) q[2];
sx q[2];
rz(1.2230011) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3712284) q[1];
sx q[1];
rz(-0.1703573) q[1];
sx q[1];
rz(-0.78070663) q[1];
rz(-pi) q[2];
rz(-0.1518941) q[3];
sx q[3];
rz(-0.23670247) q[3];
sx q[3];
rz(-3.0083002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47444433) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.8623964) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(2.352879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0448488) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(-3.0268273) q[0];
x q[1];
rz(1.6025701) q[2];
sx q[2];
rz(-1.8268158) q[2];
sx q[2];
rz(0.94871828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84856725) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(2.5679563) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2867145) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-2.9555087) q[2];
rz(0.6535334) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(2.4198467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16746178) q[0];
sx q[0];
rz(-0.92739096) q[0];
sx q[0];
rz(1.2533623) q[0];
rz(-pi) q[1];
rz(-0.53775215) q[2];
sx q[2];
rz(-2.3719412) q[2];
sx q[2];
rz(0.69307454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4221146) q[1];
sx q[1];
rz(-2.0146703) q[1];
sx q[1];
rz(1.5281954) q[1];
x q[2];
rz(-0.94354043) q[3];
sx q[3];
rz(-1.846608) q[3];
sx q[3];
rz(2.4081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.58406126) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(3.0990565) q[0];
rz(-0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-2.3775878) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873136) q[0];
sx q[0];
rz(-2.2494082) q[0];
sx q[0];
rz(-0.016088386) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6971223) q[2];
sx q[2];
rz(-0.9300803) q[2];
sx q[2];
rz(1.7510406) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9909775) q[1];
sx q[1];
rz(-2.5533278) q[1];
sx q[1];
rz(0.88612835) q[1];
rz(-2.0352827) q[3];
sx q[3];
rz(-2.1811857) q[3];
sx q[3];
rz(2.8697517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8445231) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(-2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(2.7382543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82542244) q[0];
sx q[0];
rz(-0.24992019) q[0];
sx q[0];
rz(-0.73096801) q[0];
x q[1];
rz(0.43535797) q[2];
sx q[2];
rz(-0.73409664) q[2];
sx q[2];
rz(-1.3370607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12521872) q[1];
sx q[1];
rz(-2.5902936) q[1];
sx q[1];
rz(-0.2972879) q[1];
x q[2];
rz(0.67540695) q[3];
sx q[3];
rz(-2.0007779) q[3];
sx q[3];
rz(-0.59795415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(3.1351556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2086439) q[0];
sx q[0];
rz(-0.83983487) q[0];
sx q[0];
rz(-1.9122002) q[0];
x q[1];
rz(1.8338649) q[2];
sx q[2];
rz(-0.38882133) q[2];
sx q[2];
rz(3.1281272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7240119) q[1];
sx q[1];
rz(-1.15861) q[1];
sx q[1];
rz(-0.25556232) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84900093) q[3];
sx q[3];
rz(-2.1863722) q[3];
sx q[3];
rz(-1.4403696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.6830106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4538649) q[0];
sx q[0];
rz(-1.223432) q[0];
sx q[0];
rz(-2.0485282) q[0];
x q[1];
rz(0.84476446) q[2];
sx q[2];
rz(-2.4869707) q[2];
sx q[2];
rz(2.1824333) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9142368) q[1];
sx q[1];
rz(-1.0096692) q[1];
sx q[1];
rz(-0.44740541) q[1];
rz(-1.1858995) q[3];
sx q[3];
rz(-1.2166096) q[3];
sx q[3];
rz(2.3779496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-2.356142) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.2952341) q[3];
sx q[3];
rz(2.7261962) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(-2.725504) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04836719) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(-2.918539) q[0];
rz(-pi) q[1];
rz(1.2942737) q[2];
sx q[2];
rz(-2.6575436) q[2];
sx q[2];
rz(-0.36177847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.63562288) q[1];
sx q[1];
rz(-1.9700248) q[1];
sx q[1];
rz(1.3969621) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(1.449031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-0.32315928) q[2];
rz(0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-0.51913613) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0275092) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(-1.2252349) q[0];
x q[1];
rz(-2.3245415) q[2];
sx q[2];
rz(-2.8231986) q[2];
sx q[2];
rz(1.5621834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12431006) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(-3.0148274) q[1];
rz(2.5095652) q[3];
sx q[3];
rz(-0.29187376) q[3];
sx q[3];
rz(2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.39559078) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(-1.5270773) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.4987) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258543) q[0];
sx q[0];
rz(-3.0634974) q[0];
sx q[0];
rz(1.3094835) q[0];
rz(-1.2730359) q[2];
sx q[2];
rz(-2.9794663) q[2];
sx q[2];
rz(-2.8043384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(1.0971607) q[1];
rz(0.54159553) q[3];
sx q[3];
rz(-2.9062727) q[3];
sx q[3];
rz(0.02863392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.75858086) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(2.995058) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(1.8101495) q[2];
sx q[2];
rz(-1.6806921) q[2];
sx q[2];
rz(0.42721911) q[2];
rz(0.90445789) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
