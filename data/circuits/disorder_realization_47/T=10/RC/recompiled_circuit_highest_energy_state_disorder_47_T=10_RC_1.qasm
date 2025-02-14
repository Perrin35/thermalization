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
rz(-1.9637928) q[1];
sx q[1];
rz(4.7363321) q[1];
sx q[1];
rz(7.9793032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5531122) q[0];
sx q[0];
rz(-0.40121005) q[0];
sx q[0];
rz(3.0305639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4769449) q[2];
sx q[2];
rz(-0.87021135) q[2];
sx q[2];
rz(2.075891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18516416) q[1];
sx q[1];
rz(-1.2425636) q[1];
sx q[1];
rz(0.64394711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8347391) q[3];
sx q[3];
rz(-0.05159353) q[3];
sx q[3];
rz(1.8966228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4989) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(-1.6373681) q[2];
rz(2.4046992) q[3];
sx q[3];
rz(-1.5259909) q[3];
sx q[3];
rz(-0.98114291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7351643) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(-2.6306187) q[0];
rz(-1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-2.8151292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0827507) q[0];
sx q[0];
rz(-1.3044169) q[0];
sx q[0];
rz(-0.27348117) q[0];
rz(2.9316177) q[2];
sx q[2];
rz(-2.2433503) q[2];
sx q[2];
rz(-2.5598516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.49681155) q[1];
sx q[1];
rz(-2.4678951) q[1];
sx q[1];
rz(-1.0717057) q[1];
rz(-pi) q[2];
rz(0.28387733) q[3];
sx q[3];
rz(-1.9440158) q[3];
sx q[3];
rz(2.9220327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0728077) q[0];
sx q[0];
rz(-2.7588221) q[0];
sx q[0];
rz(-2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-0.41789564) q[1];
sx q[1];
rz(0.03104041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22716534) q[0];
sx q[0];
rz(-1.7591159) q[0];
sx q[0];
rz(0.31115599) q[0];
x q[1];
rz(0.91340567) q[2];
sx q[2];
rz(-1.3904461) q[2];
sx q[2];
rz(-3.0424398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3794897) q[1];
sx q[1];
rz(-0.65564686) q[1];
sx q[1];
rz(-1.9807705) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1150041) q[3];
sx q[3];
rz(-2.5864161) q[3];
sx q[3];
rz(2.7257277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62640181) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(-0.82702965) q[2];
rz(0.51860297) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(-1.5145068) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14722918) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(1.9299141) q[0];
rz(1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(-1.3406219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692684) q[0];
sx q[0];
rz(-2.4256673) q[0];
sx q[0];
rz(0.47874079) q[0];
x q[1];
rz(-0.2366613) q[2];
sx q[2];
rz(-1.4259031) q[2];
sx q[2];
rz(-1.7828816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0405492) q[1];
sx q[1];
rz(-1.2770997) q[1];
sx q[1];
rz(-2.5459176) q[1];
rz(2.2068437) q[3];
sx q[3];
rz(-2.6451561) q[3];
sx q[3];
rz(-1.8938586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(-1.1411257) q[2];
rz(1.6107669) q[3];
sx q[3];
rz(-0.96735668) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56213266) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(0.60254565) q[0];
rz(-2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(2.3604732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8373708) q[0];
sx q[0];
rz(-1.8580372) q[0];
sx q[0];
rz(2.9648613) q[0];
rz(-0.94046142) q[2];
sx q[2];
rz(-0.9113208) q[2];
sx q[2];
rz(2.3465921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32051099) q[1];
sx q[1];
rz(-0.85144224) q[1];
sx q[1];
rz(3.0266552) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72569287) q[3];
sx q[3];
rz(-1.7039434) q[3];
sx q[3];
rz(0.6729047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96131229) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(-2.9126634) q[2];
rz(-1.3198352) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22885403) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(1.1786906) q[0];
rz(-1.9121869) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(-1.8454525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4018009) q[0];
sx q[0];
rz(-1.290088) q[0];
sx q[0];
rz(-0.14174353) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4086094) q[2];
sx q[2];
rz(-2.146581) q[2];
sx q[2];
rz(-2.1819262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3553222) q[1];
sx q[1];
rz(-2.4679524) q[1];
sx q[1];
rz(-2.4060529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9518106) q[3];
sx q[3];
rz(-2.6468607) q[3];
sx q[3];
rz(-3.1056946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.487315) q[2];
sx q[2];
rz(-2.2389905) q[2];
sx q[2];
rz(-0.26228341) q[2];
rz(0.30019635) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009785) q[0];
sx q[0];
rz(-2.746026) q[0];
sx q[0];
rz(1.8390919) q[0];
rz(-0.66849661) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(2.1065333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9205017) q[0];
sx q[0];
rz(-1.1701487) q[0];
sx q[0];
rz(2.9095838) q[0];
x q[1];
rz(-2.9547353) q[2];
sx q[2];
rz(-1.2070939) q[2];
sx q[2];
rz(-1.1101369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25495249) q[1];
sx q[1];
rz(-1.4952085) q[1];
sx q[1];
rz(1.2264538) q[1];
rz(-pi) q[2];
rz(-0.74372282) q[3];
sx q[3];
rz(-1.111314) q[3];
sx q[3];
rz(2.2032025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11789007) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(1.6560076) q[2];
rz(-2.8591136) q[3];
sx q[3];
rz(-1.7829203) q[3];
sx q[3];
rz(-0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(-1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(-2.5801632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0208541) q[0];
sx q[0];
rz(-2.0971813) q[0];
sx q[0];
rz(0.83197439) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8738633) q[2];
sx q[2];
rz(-0.91721877) q[2];
sx q[2];
rz(-0.67959626) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16374396) q[1];
sx q[1];
rz(-1.9789961) q[1];
sx q[1];
rz(-0.331649) q[1];
rz(1.9188324) q[3];
sx q[3];
rz(-1.819918) q[3];
sx q[3];
rz(3.0580229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5795035) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(1.7378463) q[2];
rz(1.7956519) q[3];
sx q[3];
rz(-2.3965049) q[3];
sx q[3];
rz(-2.6534206) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3812934) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(-2.2316933) q[0];
rz(-1.1970041) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(2.2534456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040374856) q[0];
sx q[0];
rz(-1.4618317) q[0];
sx q[0];
rz(0.094281406) q[0];
x q[1];
rz(1.6309225) q[2];
sx q[2];
rz(-0.91183582) q[2];
sx q[2];
rz(1.9163641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8348114) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(2.4624834) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1040583) q[3];
sx q[3];
rz(-1.2753295) q[3];
sx q[3];
rz(-2.7151544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1278648) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(-2.0590651) q[2];
rz(2.9108289) q[3];
sx q[3];
rz(-1.328238) q[3];
sx q[3];
rz(-2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.3103264) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(2.5600774) q[0];
rz(2.5260018) q[1];
sx q[1];
rz(-0.94771996) q[1];
sx q[1];
rz(-1.8448578) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0688419) q[0];
sx q[0];
rz(-1.5525991) q[0];
sx q[0];
rz(-0.010404603) q[0];
rz(1.5611951) q[2];
sx q[2];
rz(-1.9148972) q[2];
sx q[2];
rz(0.32479686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2263068) q[1];
sx q[1];
rz(-1.9939594) q[1];
sx q[1];
rz(-2.8598815) q[1];
x q[2];
rz(2.6526645) q[3];
sx q[3];
rz(-2.4081514) q[3];
sx q[3];
rz(-1.5239925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41946188) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(0.44931832) q[2];
rz(-2.8214473) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(1.2464574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768059) q[0];
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
rz(2.7560227) q[3];
sx q[3];
rz(-1.1299993) q[3];
sx q[3];
rz(-0.53408505) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
