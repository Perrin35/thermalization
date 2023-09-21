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
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0770744) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(1.6032739) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50813714) q[2];
sx q[2];
rz(-1.7388441) q[2];
sx q[2];
rz(0.44067581) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5832311) q[1];
sx q[1];
rz(-1.6915295) q[1];
sx q[1];
rz(1.4503149) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5343127) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46368018) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(0.15788831) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(0.78871361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4855027) q[0];
sx q[0];
rz(-0.1682818) q[0];
sx q[0];
rz(-2.3165354) q[0];
rz(-1.6025701) q[2];
sx q[2];
rz(-1.3147768) q[2];
sx q[2];
rz(-2.1928744) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2930254) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(-2.5679563) q[1];
rz(-pi) q[2];
rz(2.2867145) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(2.9555087) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(2.511456) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-2.4198467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8080224) q[0];
sx q[0];
rz(-2.434242) q[0];
sx q[0];
rz(0.39444123) q[0];
rz(-pi) q[1];
rz(-2.0314991) q[2];
sx q[2];
rz(-0.93020541) q[2];
sx q[2];
rz(5.0355807e-05) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5210515) q[1];
sx q[1];
rz(-2.6958145) q[1];
sx q[1];
rz(3.0522703) q[1];
rz(-pi) q[2];
rz(2.8053022) q[3];
sx q[3];
rz(-2.170917) q[3];
sx q[3];
rz(-2.1093413) q[3];
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
rz(-2.2325113) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
sx q[1];
rz(-0.50023729) q[1];
sx q[1];
rz(0.76400486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9616868) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(1.5907445) q[0];
rz(2.4970384) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(-2.8855756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92661392) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(-2.742393) q[1];
rz(-pi) q[2];
x q[2];
rz(1.10631) q[3];
sx q[3];
rz(-0.96040695) q[3];
sx q[3];
rz(-2.8697517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9168636) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(-2.9149756) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(0.40333834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628495) q[0];
sx q[0];
rz(-1.3855977) q[0];
sx q[0];
rz(-1.7395822) q[0];
rz(-1.9344159) q[2];
sx q[2];
rz(-0.91797963) q[2];
sx q[2];
rz(-1.8967241) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4407318) q[1];
sx q[1];
rz(-1.4167538) q[1];
sx q[1];
rz(-2.6101019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5087318) q[3];
sx q[3];
rz(-2.359458) q[3];
sx q[3];
rz(2.6485505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(0.0064370357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9329487) q[0];
sx q[0];
rz(-0.83983487) q[0];
sx q[0];
rz(-1.9122002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1941031) q[2];
sx q[2];
rz(-1.4720535) q[2];
sx q[2];
rz(-1.313098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2575063) q[1];
sx q[1];
rz(-1.8045366) q[1];
sx q[1];
rz(-1.9952378) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75140679) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-0.71099647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(-0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(0.78272351) q[1];
sx q[1];
rz(-2.0109773) q[1];
sx q[1];
rz(1.9810716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4538649) q[0];
sx q[0];
rz(-1.223432) q[0];
sx q[0];
rz(2.0485282) q[0];
x q[1];
rz(1.0497401) q[2];
sx q[2];
rz(-1.9869291) q[2];
sx q[2];
rz(-0.001948826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50642636) q[1];
sx q[1];
rz(-2.4392358) q[1];
sx q[1];
rz(-0.96794767) q[1];
rz(0.79302391) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(-1.6263863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(2.356142) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(1.998741) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4131665) q[0];
sx q[0];
rz(-2.6217555) q[0];
sx q[0];
rz(1.1632989) q[0];
rz(-pi) q[1];
rz(1.1025238) q[2];
sx q[2];
rz(-1.6981914) q[2];
sx q[2];
rz(0.96291908) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0033274) q[1];
sx q[1];
rz(-1.4107553) q[1];
sx q[1];
rz(0.40469594) q[1];
rz(-pi) q[2];
rz(-0.11085005) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(-1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828744) q[0];
sx q[0];
rz(-1.9472194) q[0];
sx q[0];
rz(2.9988078) q[0];
x q[1];
rz(0.81705117) q[2];
sx q[2];
rz(-2.8231986) q[2];
sx q[2];
rz(1.5621834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6053033) q[1];
sx q[1];
rz(-1.481206) q[1];
sx q[1];
rz(-2.3593966) q[1];
x q[2];
rz(2.5095652) q[3];
sx q[3];
rz(-2.8497189) q[3];
sx q[3];
rz(0.29575086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460019) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.4987) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1536627) q[0];
sx q[0];
rz(-1.6462353) q[0];
sx q[0];
rz(-3.1213785) q[0];
rz(-pi) q[1];
rz(3.0936436) q[2];
sx q[2];
rz(-1.4158632) q[2];
sx q[2];
rz(-0.63873728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1415256) q[1];
sx q[1];
rz(-0.86874092) q[1];
sx q[1];
rz(-2.6932004) q[1];
rz(-pi) q[2];
rz(0.20262952) q[3];
sx q[3];
rz(-1.6912795) q[3];
sx q[3];
rz(2.1287624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.8101495) q[2];
sx q[2];
rz(-1.4609006) q[2];
sx q[2];
rz(-2.7143735) q[2];
rz(-2.2371348) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
