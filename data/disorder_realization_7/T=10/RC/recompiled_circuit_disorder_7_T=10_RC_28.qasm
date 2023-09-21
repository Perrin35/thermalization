OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1768567) q[0];
sx q[0];
rz(-0.86369792) q[0];
sx q[0];
rz(0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62030828) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(1.3814397) q[0];
rz(2.7896499) q[2];
sx q[2];
rz(-1.8004187) q[2];
sx q[2];
rz(0.72598347) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37571733) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(2.767763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8455891) q[3];
sx q[3];
rz(-1.4010251) q[3];
sx q[3];
rz(0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-0.10736297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023022368) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(0.9704216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41610833) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(-2.1925418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0432942) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(2.1143101) q[1];
rz(-0.3932088) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(-1.8821017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28856453) q[0];
sx q[0];
rz(-1.0596501) q[0];
sx q[0];
rz(1.9938064) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9627377) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(-2.4286963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(2.9134681) q[1];
rz(-pi) q[2];
rz(2.9186967) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(-0.086829348) q[0];
rz(2.343802) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(2.4385902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9490337) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(2.2662524) q[1];
x q[2];
rz(0.66588464) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(-0.22660412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.8866084) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(0.19651861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5614583) q[0];
sx q[0];
rz(-2.723454) q[0];
sx q[0];
rz(2.3925376) q[0];
rz(-pi) q[1];
rz(-0.032012149) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(0.33547685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4938441) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(2.8791048) q[1];
x q[2];
rz(-1.7408095) q[3];
sx q[3];
rz(-1.2983054) q[3];
sx q[3];
rz(-1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(-2.939558) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.3938168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73035875) q[0];
sx q[0];
rz(-2.8011836) q[0];
sx q[0];
rz(2.911724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9453085) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(1.166677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-2.1761314) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4059541) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(-2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-0.21025118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41482571) q[0];
sx q[0];
rz(-1.6370156) q[0];
sx q[0];
rz(1.7050752) q[0];
rz(0.37961752) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47760689) q[1];
sx q[1];
rz(-0.37316445) q[1];
sx q[1];
rz(0.54877703) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(-1.9853026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.7604527) q[2];
rz(0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64637) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(-0.97672021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57709549) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(-0.45678267) q[0];
rz(-pi) q[1];
x q[1];
rz(3.061053) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(1.5921519) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50423065) q[1];
sx q[1];
rz(-1.1334051) q[1];
sx q[1];
rz(-2.7573542) q[1];
x q[2];
rz(0.85985698) q[3];
sx q[3];
rz(-2.007774) q[3];
sx q[3];
rz(-2.6881998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-0.92064944) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9793648) q[0];
sx q[0];
rz(-0.75408903) q[0];
sx q[0];
rz(0.92569949) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1295206) q[2];
sx q[2];
rz(-1.8555102) q[2];
sx q[2];
rz(0.65288359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2194849) q[1];
sx q[1];
rz(-1.9003632) q[1];
sx q[1];
rz(0.12164128) q[1];
rz(0.00021342834) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(0.40094024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(0.10458874) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(1.1178281) q[0];
rz(-2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871646) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(0.54425311) q[0];
x q[1];
rz(1.2043578) q[2];
sx q[2];
rz(-1.2841184) q[2];
sx q[2];
rz(-0.8358801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6798903) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(0.093613503) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0685365) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(1.6554228) q[2];
sx q[2];
rz(-1.6783236) q[2];
sx q[2];
rz(2.8477737) q[2];
rz(1.312064) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
