OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(-0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(2.2367509) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358937) q[0];
sx q[0];
rz(-2.4568757) q[0];
sx q[0];
rz(-0.84795714) q[0];
rz(-pi) q[1];
rz(-2.4279847) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(-1.7507391) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(0.25524615) q[1];
x q[2];
rz(-0.18890394) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(-0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(3.0453483) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-2.3764215) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-2.4786425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11766079) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(-0.06225417) q[0];
rz(-pi) q[1];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1124277) q[1];
sx q[1];
rz(-1.4883092) q[1];
sx q[1];
rz(-2.6998181) q[1];
rz(-pi) q[2];
rz(-2.8421721) q[3];
sx q[3];
rz(-2.6357108) q[3];
sx q[3];
rz(-0.29380709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(0.66550231) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(-1.3177419) q[3];
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
rz(pi/2) q[0];
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
rz(-0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(0.51868784) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.737239) q[0];
sx q[0];
rz(-0.81439942) q[0];
sx q[0];
rz(1.0712207) q[0];
rz(-pi) q[1];
x q[1];
rz(1.475004) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(-1.4694627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7171214) q[1];
sx q[1];
rz(-0.79454225) q[1];
sx q[1];
rz(-0.24309991) q[1];
rz(-pi) q[2];
rz(-0.014702602) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1888694) q[0];
sx q[0];
rz(-2.1892622) q[0];
sx q[0];
rz(-1.6678715) q[0];
x q[1];
rz(2.3344343) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(-1.6023139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0780371) q[1];
sx q[1];
rz(-2.6647898) q[1];
sx q[1];
rz(0.30492353) q[1];
rz(-1.7771878) q[3];
sx q[3];
rz(-1.5406113) q[3];
sx q[3];
rz(-0.15011945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84232932) q[0];
sx q[0];
rz(-3.0780601) q[0];
sx q[0];
rz(-2.1701943) q[0];
rz(2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(-2.9074557) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7580326) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(0.010096117) q[1];
rz(1.1957937) q[3];
sx q[3];
rz(-0.71205322) q[3];
sx q[3];
rz(-0.28111162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(-2.8862254) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(-0.52608144) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(0.88476673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381829) q[0];
sx q[0];
rz(-0.38823715) q[0];
sx q[0];
rz(1.7675179) q[0];
rz(-pi) q[1];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5846328) q[2];
sx q[2];
rz(-2.7163598) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92237597) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(-2.8273696) q[1];
rz(-pi) q[2];
rz(-0.24598908) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46985627) q[0];
sx q[0];
rz(-2.0977019) q[0];
sx q[0];
rz(-2.0335474) q[0];
x q[1];
rz(1.0211208) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(2.6598425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43436189) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-1.93768) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(2.4550408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(0.7094267) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(0.27871305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3942791) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-0.5745116) q[0];
rz(-pi) q[1];
rz(-1.3271354) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.3512163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7396001) q[1];
sx q[1];
rz(-2.7584689) q[1];
sx q[1];
rz(-0.95462228) q[1];
rz(-pi) q[2];
rz(1.2392427) q[3];
sx q[3];
rz(-0.23570508) q[3];
sx q[3];
rz(-1.9086259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(-2.6682207) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.8737428) q[0];
x q[1];
rz(1.5683453) q[2];
sx q[2];
rz(-0.8674538) q[2];
sx q[2];
rz(-3.0424812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0513623) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(2.5999703) q[1];
rz(-pi) q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(0.70860329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(2.4401869) q[0];
rz(-2.2258863) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(-1.2385626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68574821) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(2.8932163) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5194703) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(-2.3527956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.726767) q[1];
sx q[1];
rz(-1.4531724) q[1];
sx q[1];
rz(-2.2284501) q[1];
rz(-1.402461) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(-1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713521) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(-1.3958037) q[3];
sx q[3];
rz(-0.77944118) q[3];
sx q[3];
rz(-0.13164095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
