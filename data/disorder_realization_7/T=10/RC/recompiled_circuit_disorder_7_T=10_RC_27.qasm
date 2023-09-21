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
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(-0.62682682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1683567) q[0];
sx q[0];
rz(-1.3827948) q[0];
sx q[0];
rz(-0.12113916) q[0];
rz(-pi) q[1];
rz(-0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(-0.28996224) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1605692) q[1];
sx q[1];
rz(-1.261928) q[1];
sx q[1];
rz(0.94998756) q[1];
x q[2];
rz(-1.2960035) q[3];
sx q[3];
rz(-1.4010251) q[3];
sx q[3];
rz(0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(0.75479341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24793808) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0221314) q[0];
sx q[0];
rz(-2.4298926) q[0];
sx q[0];
rz(0.91711451) q[0];
rz(0.58789247) q[2];
sx q[2];
rz(-1.3304454) q[2];
sx q[2];
rz(-2.8628778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0432942) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(-2.1143101) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7530305) q[3];
sx q[3];
rz(-1.4149168) q[3];
sx q[3];
rz(0.673783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(-2.4798933) q[2];
rz(0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(2.0592164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6425991) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(-0.55143349) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9627377) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(-0.71289635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2228959) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(3.0547633) q[0];
rz(-2.5627665) q[2];
sx q[2];
rz(-2.0818713) q[2];
sx q[2];
rz(1.5103112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9490337) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(2.2662524) q[1];
rz(-0.66588464) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(-2.9149885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-0.19651861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949961) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(-2.8269935) q[0];
rz(-pi) q[1];
rz(2.0318803) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(-1.2210786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4938441) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(-0.26248787) q[1];
rz(2.8653141) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(-2.7460263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(-0.20203461) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(-1.3938168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73035875) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(-2.911724) q[0];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-0.44484777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5382232) q[1];
sx q[1];
rz(-1.9586316) q[1];
sx q[1];
rz(2.8549854) q[1];
x q[2];
rz(-0.73563852) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(-2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(-1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(0.21025118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41482571) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(1.4365175) q[0];
rz(-pi) q[1];
rz(0.37961752) q[2];
sx q[2];
rz(-0.175975) q[2];
sx q[2];
rz(1.1773674) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10340313) q[1];
sx q[1];
rz(-1.8870755) q[1];
sx q[1];
rz(1.7722539) q[1];
rz(-0.85150163) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(-2.7281318) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(2.1648724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8551089) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(2.8141862) q[0];
rz(1.7325399) q[2];
sx q[2];
rz(-1.6502893) q[2];
sx q[2];
rz(0.0083991945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2663914) q[1];
sx q[1];
rz(-2.5677263) q[1];
sx q[1];
rz(0.89504524) q[1];
rz(-2.1920491) q[3];
sx q[3];
rz(-0.81406677) q[3];
sx q[3];
rz(0.66063389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(0.92064944) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-2.2424973) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16222787) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(2.2158932) q[0];
rz(-1.012072) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(-0.65288359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60914674) q[1];
sx q[1];
rz(-1.6858613) q[1];
sx q[1];
rz(1.902641) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5713463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(0.88275638) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1966232) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(2.1150132) q[0];
rz(-0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(-3.0416833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9576479) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(-1.8326879) q[1];
x q[2];
rz(-2.2393353) q[3];
sx q[3];
rz(-0.60562953) q[3];
sx q[3];
rz(2.827364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(1.2697521) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-0.10791049) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
rz(3.0222223) q[3];
sx q[3];
rz(-1.8277677) q[3];
sx q[3];
rz(1.2895332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];