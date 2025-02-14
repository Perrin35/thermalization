OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90836877) q[0];
sx q[0];
rz(-1.5923192) q[0];
sx q[0];
rz(0.17750658) q[0];
rz(-1.5818051) q[1];
sx q[1];
rz(-1.3270562) q[1];
sx q[1];
rz(1.1787193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27025199) q[0];
sx q[0];
rz(-2.1553023) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(2.841973) q[2];
sx q[2];
rz(-2.1914406) q[2];
sx q[2];
rz(-2.621068) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.385966) q[1];
sx q[1];
rz(-2.3896791) q[1];
sx q[1];
rz(-2.2346943) q[1];
rz(-pi) q[2];
rz(0.13362515) q[3];
sx q[3];
rz(-0.84098626) q[3];
sx q[3];
rz(-1.3624024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0796279) q[2];
sx q[2];
rz(-2.2140333) q[2];
sx q[2];
rz(2.8947042) q[2];
rz(0.43761349) q[3];
sx q[3];
rz(-2.103431) q[3];
sx q[3];
rz(-0.26052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058411773) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(2.0368982) q[0];
rz(-0.80766922) q[1];
sx q[1];
rz(-0.75391114) q[1];
sx q[1];
rz(0.63899904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7154626) q[0];
sx q[0];
rz(-1.3952266) q[0];
sx q[0];
rz(2.8605754) q[0];
rz(0.27490567) q[2];
sx q[2];
rz(-1.5228809) q[2];
sx q[2];
rz(2.3570182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39589992) q[1];
sx q[1];
rz(-1.2077161) q[1];
sx q[1];
rz(-2.3895742) q[1];
rz(0.86231972) q[3];
sx q[3];
rz(-0.42565027) q[3];
sx q[3];
rz(-1.6274393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4421926) q[2];
sx q[2];
rz(-1.7510119) q[2];
sx q[2];
rz(-2.3244582) q[2];
rz(-2.6160431) q[3];
sx q[3];
rz(-2.3216129) q[3];
sx q[3];
rz(1.1901963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84080559) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(-2.9817885) q[0];
rz(-2.6218759) q[1];
sx q[1];
rz(-1.0451885) q[1];
sx q[1];
rz(-2.2249178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063112886) q[0];
sx q[0];
rz(-1.0005299) q[0];
sx q[0];
rz(-1.6018014) q[0];
rz(1.3266852) q[2];
sx q[2];
rz(-2.296431) q[2];
sx q[2];
rz(0.47032088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0708051) q[1];
sx q[1];
rz(-0.19529058) q[1];
sx q[1];
rz(1.8205804) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5289747) q[3];
sx q[3];
rz(-1.4605986) q[3];
sx q[3];
rz(-1.8511021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89797574) q[2];
sx q[2];
rz(-2.6349082) q[2];
sx q[2];
rz(-1.766073) q[2];
rz(-0.033179387) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(0.83056617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0783949) q[0];
sx q[0];
rz(-2.4802408) q[0];
sx q[0];
rz(0.98264328) q[0];
rz(2.4962418) q[1];
sx q[1];
rz(-1.9159562) q[1];
sx q[1];
rz(2.7243848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010734) q[0];
sx q[0];
rz(-0.60581453) q[0];
sx q[0];
rz(-2.5941501) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9511767) q[2];
sx q[2];
rz(-1.6350758) q[2];
sx q[2];
rz(-0.96997875) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4667838) q[1];
sx q[1];
rz(-2.4655113) q[1];
sx q[1];
rz(-3.0463808) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44261114) q[3];
sx q[3];
rz(-2.5472982) q[3];
sx q[3];
rz(-1.3292918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30503201) q[2];
sx q[2];
rz(-2.8053668) q[2];
sx q[2];
rz(-1.1870556) q[2];
rz(0.120397) q[3];
sx q[3];
rz(-1.4056561) q[3];
sx q[3];
rz(2.9482237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(2.8315941) q[0];
sx q[0];
rz(-3.031142) q[0];
sx q[0];
rz(2.4647392) q[0];
rz(1.5325158) q[1];
sx q[1];
rz(-2.0307505) q[1];
sx q[1];
rz(-2.9422876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4368962) q[0];
sx q[0];
rz(-1.4802209) q[0];
sx q[0];
rz(1.3170088) q[0];
x q[1];
rz(-1.6930313) q[2];
sx q[2];
rz(-1.2265556) q[2];
sx q[2];
rz(-0.93512541) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90016323) q[1];
sx q[1];
rz(-1.8668264) q[1];
sx q[1];
rz(-2.7524292) q[1];
rz(1.5478743) q[3];
sx q[3];
rz(-2.8646152) q[3];
sx q[3];
rz(-2.2131526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4946332) q[2];
sx q[2];
rz(-0.55611742) q[2];
sx q[2];
rz(2.6456918) q[2];
rz(-0.9209218) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(0.31057772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50716901) q[0];
sx q[0];
rz(-1.2260219) q[0];
sx q[0];
rz(0.41967151) q[0];
rz(-0.93027973) q[1];
sx q[1];
rz(-2.1280961) q[1];
sx q[1];
rz(-1.0118265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9032325) q[0];
sx q[0];
rz(-1.6201265) q[0];
sx q[0];
rz(-0.65261638) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18380766) q[2];
sx q[2];
rz(-1.70924) q[2];
sx q[2];
rz(-0.02414298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25988425) q[1];
sx q[1];
rz(-1.0233423) q[1];
sx q[1];
rz(-0.24357067) q[1];
rz(-1.3023754) q[3];
sx q[3];
rz(-0.61385775) q[3];
sx q[3];
rz(0.74210244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73152995) q[2];
sx q[2];
rz(-2.5750934) q[2];
sx q[2];
rz(-2.0274963) q[2];
rz(-1.3725494) q[3];
sx q[3];
rz(-1.0053582) q[3];
sx q[3];
rz(0.70926386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0396742) q[0];
sx q[0];
rz(-1.8721975) q[0];
sx q[0];
rz(1.6495548) q[0];
rz(0.44293013) q[1];
sx q[1];
rz(-1.8699173) q[1];
sx q[1];
rz(-2.2992004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5076601) q[0];
sx q[0];
rz(-2.2325987) q[0];
sx q[0];
rz(-2.4217601) q[0];
x q[1];
rz(0.90656735) q[2];
sx q[2];
rz(-1.7822872) q[2];
sx q[2];
rz(1.6935401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2017224) q[1];
sx q[1];
rz(-2.2186167) q[1];
sx q[1];
rz(0.52380162) q[1];
rz(-0.042875127) q[3];
sx q[3];
rz(-1.6742286) q[3];
sx q[3];
rz(-2.7989911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3288021) q[2];
sx q[2];
rz(-2.2219358) q[2];
sx q[2];
rz(0.021050464) q[2];
rz(0.99714315) q[3];
sx q[3];
rz(-0.54233426) q[3];
sx q[3];
rz(1.5525345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7954623) q[0];
sx q[0];
rz(-1.902245) q[0];
sx q[0];
rz(-1.963266) q[0];
rz(1.0109674) q[1];
sx q[1];
rz(-1.6594454) q[1];
sx q[1];
rz(-2.7706026) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8218433) q[0];
sx q[0];
rz(-0.25938636) q[0];
sx q[0];
rz(0.46617561) q[0];
rz(-1.402032) q[2];
sx q[2];
rz(-1.4764364) q[2];
sx q[2];
rz(0.33359087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4996262) q[1];
sx q[1];
rz(-1.2807944) q[1];
sx q[1];
rz(-2.9846342) q[1];
rz(-pi) q[2];
rz(0.31412394) q[3];
sx q[3];
rz(-1.1289136) q[3];
sx q[3];
rz(-1.266138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79169881) q[2];
sx q[2];
rz(-2.0588304) q[2];
sx q[2];
rz(0.91378158) q[2];
rz(2.6036116) q[3];
sx q[3];
rz(-2.3086083) q[3];
sx q[3];
rz(-0.37916455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.051788483) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(2.9096933) q[0];
rz(1.505835) q[1];
sx q[1];
rz(-1.4488522) q[1];
sx q[1];
rz(2.4969782) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5568791) q[0];
sx q[0];
rz(-1.9125835) q[0];
sx q[0];
rz(0.36447592) q[0];
x q[1];
rz(1.5455075) q[2];
sx q[2];
rz(-2.1754258) q[2];
sx q[2];
rz(2.189332) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48460173) q[1];
sx q[1];
rz(-1.4047523) q[1];
sx q[1];
rz(-3.1216122) q[1];
x q[2];
rz(0.70213153) q[3];
sx q[3];
rz(-1.9117179) q[3];
sx q[3];
rz(0.11521712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5171648) q[2];
sx q[2];
rz(-1.6021873) q[2];
sx q[2];
rz(1.2702764) q[2];
rz(0.43954784) q[3];
sx q[3];
rz(-0.63616532) q[3];
sx q[3];
rz(-2.6308681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9257833) q[0];
sx q[0];
rz(-1.2405115) q[0];
sx q[0];
rz(-1.0260169) q[0];
rz(2.6846474) q[1];
sx q[1];
rz(-0.89070717) q[1];
sx q[1];
rz(-2.2073726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608604) q[0];
sx q[0];
rz(-2.6965158) q[0];
sx q[0];
rz(-2.829192) q[0];
rz(-pi) q[1];
x q[1];
rz(3.031557) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(2.2243654) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10623744) q[1];
sx q[1];
rz(-0.83939531) q[1];
sx q[1];
rz(1.3421571) q[1];
rz(-2.2226376) q[3];
sx q[3];
rz(-2.0522293) q[3];
sx q[3];
rz(1.0713861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.702174) q[2];
sx q[2];
rz(-2.4472523) q[2];
sx q[2];
rz(1.0609421) q[2];
rz(0.98961467) q[3];
sx q[3];
rz(-1.7695534) q[3];
sx q[3];
rz(2.8770679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2117352) q[0];
sx q[0];
rz(-2.0462357) q[0];
sx q[0];
rz(-0.76475058) q[0];
rz(-1.9542971) q[1];
sx q[1];
rz(-1.7056414) q[1];
sx q[1];
rz(2.5797226) q[1];
rz(-1.622772) q[2];
sx q[2];
rz(-1.5511654) q[2];
sx q[2];
rz(1.8775107) q[2];
rz(2.3791947) q[3];
sx q[3];
rz(-2.5350606) q[3];
sx q[3];
rz(1.4192941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
