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
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(-1.5919332) q[1];
sx q[1];
rz(-3.0631493) q[1];
sx q[1];
rz(0.6426386) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064518236) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(1.5383188) q[0];
x q[1];
rz(0.3354934) q[2];
sx q[2];
rz(-2.6087084) q[2];
sx q[2];
rz(0.83836183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77036422) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(2.9075165) q[3];
sx q[3];
rz(-1.6062859) q[3];
sx q[3];
rz(-1.8518098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.2791963) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46368018) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-0.78871361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4855027) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(-2.3165354) q[0];
rz(-pi) q[1];
rz(0.2561432) q[2];
sx q[2];
rz(-1.5400585) q[2];
sx q[2];
rz(-0.61402938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5309108) q[1];
sx q[1];
rz(-2.1149939) q[1];
sx q[1];
rz(-1.2128085) q[1];
rz(-0.7699645) q[3];
sx q[3];
rz(-1.0125481) q[3];
sx q[3];
rz(-1.6127197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.366189) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(2.4198467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5979413) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(0.6681722) q[0];
rz(-pi) q[1];
rz(-2.4475054) q[2];
sx q[2];
rz(-1.9352479) q[2];
sx q[2];
rz(1.2823766) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0112146) q[1];
sx q[1];
rz(-1.6092669) q[1];
sx q[1];
rz(2.6973666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94354043) q[3];
sx q[3];
rz(-1.846608) q[3];
sx q[3];
rz(-0.73345473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(0.90908137) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(-2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(-0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-2.3775878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1799058) q[0];
sx q[0];
rz(-2.4628203) q[0];
sx q[0];
rz(-1.5907445) q[0];
rz(1.6971223) q[2];
sx q[2];
rz(-2.2115123) q[2];
sx q[2];
rz(1.7510406) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2149787) q[1];
sx q[1];
rz(-2.0151295) q[1];
sx q[1];
rz(2.742393) q[1];
x q[2];
rz(-0.56960168) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(-0.99018712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8445231) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(2.9149756) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(2.7382543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82542244) q[0];
sx q[0];
rz(-2.8916725) q[0];
sx q[0];
rz(-0.73096801) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2071768) q[2];
sx q[2];
rz(-0.91797963) q[2];
sx q[2];
rz(1.2448685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4407318) q[1];
sx q[1];
rz(-1.4167538) q[1];
sx q[1];
rz(0.53149077) q[1];
x q[2];
rz(0.63286085) q[3];
sx q[3];
rz(-2.359458) q[3];
sx q[3];
rz(0.49304214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.430442) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(0.0064370357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.54654) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(2.3810054) q[0];
rz(-1.8338649) q[2];
sx q[2];
rz(-2.7527713) q[2];
sx q[2];
rz(-0.013465492) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16033123) q[1];
sx q[1];
rz(-2.6604974) q[1];
sx q[1];
rz(-2.0950003) q[1];
rz(-2.3859343) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(-0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7727938) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(-2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098175123) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(-3.085882) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502064) q[0];
sx q[0];
rz(-1.123748) q[0];
sx q[0];
rz(-2.7544751) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2968282) q[2];
sx q[2];
rz(-2.4869707) q[2];
sx q[2];
rz(2.1824333) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5481422) q[1];
sx q[1];
rz(-1.1957809) q[1];
sx q[1];
rz(2.1795991) q[1];
rz(-2.3485687) q[3];
sx q[3];
rz(-0.51699713) q[3];
sx q[3];
rz(-1.5152064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(2.356142) q[2];
rz(2.3857332) q[3];
sx q[3];
rz(-1.2952341) q[3];
sx q[3];
rz(-2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(-0.41608861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7284262) q[0];
sx q[0];
rz(-2.6217555) q[0];
sx q[0];
rz(1.1632989) q[0];
rz(0.14256723) q[2];
sx q[2];
rz(-1.1066184) q[2];
sx q[2];
rz(0.67205059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21129747) q[1];
sx q[1];
rz(-2.7080309) q[1];
sx q[1];
rz(-2.7525206) q[1];
x q[2];
rz(-3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(0.51913613) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3648758) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(-1.1908635) q[0];
rz(-pi) q[1];
rz(-0.22186188) q[2];
sx q[2];
rz(-1.8010745) q[2];
sx q[2];
rz(-0.78267539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.054159315) q[1];
sx q[1];
rz(-2.3490153) q[1];
sx q[1];
rz(-1.6967609) q[1];
rz(-pi) q[2];
rz(-0.63202745) q[3];
sx q[3];
rz(-2.8497189) q[3];
sx q[3];
rz(-2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-2.273902) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-0.51122558) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.4987) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98793) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(-3.1213785) q[0];
rz(1.8685568) q[2];
sx q[2];
rz(-2.9794663) q[2];
sx q[2];
rz(-2.8043384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1415256) q[1];
sx q[1];
rz(-0.86874092) q[1];
sx q[1];
rz(-2.6932004) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6937709) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(-2.995058) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-2.0064034) q[2];
sx q[2];
rz(-0.26293593) q[2];
sx q[2];
rz(1.5756366) q[2];
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
