OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456822) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(2.1366871) q[0];
rz(0.35194273) q[2];
sx q[2];
rz(-1.8004187) q[2];
sx q[2];
rz(2.4156092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1605692) q[1];
sx q[1];
rz(-1.261928) q[1];
sx q[1];
rz(0.94998756) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9653373) q[3];
sx q[3];
rz(-1.8415383) q[3];
sx q[3];
rz(0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-0.75479341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(2.1564116) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-3.0342297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1185703) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(0.9704216) q[0];
rz(-1.2843578) q[2];
sx q[2];
rz(-1.0019433) q[2];
sx q[2];
rz(1.4494277) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0432942) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(-1.0272825) q[1];
rz(1.739005) q[3];
sx q[3];
rz(-1.9544) q[3];
sx q[3];
rz(-2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-2.4798933) q[2];
rz(3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6425991) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(-2.5901592) q[0];
rz(-1.1788549) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(-0.71289635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12990002) q[1];
sx q[1];
rz(-1.3492279) q[1];
sx q[1];
rz(1.3264015) q[1];
rz(-pi) q[2];
rz(2.9186967) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(2.1427758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(0.48809537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29823829) q[0];
sx q[0];
rz(-2.6697864) q[0];
sx q[0];
rz(1.3993553) q[0];
rz(-2.5627665) q[2];
sx q[2];
rz(-2.0818713) q[2];
sx q[2];
rz(1.5103112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9490337) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(0.87534027) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5591764) q[3];
sx q[3];
rz(-1.9786668) q[3];
sx q[3];
rz(1.2553314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-2.945074) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949961) q[0];
sx q[0];
rz(-1.8509522) q[0];
sx q[0];
rz(2.8269935) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5064266) q[2];
sx q[2];
rz(-0.46191051) q[2];
sx q[2];
rz(-0.40735746) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6424751) q[1];
sx q[1];
rz(-2.827008) q[1];
sx q[1];
rz(-2.5423074) q[1];
x q[2];
rz(2.8653141) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(-2.7460263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6348811) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48702792) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(1.4902671) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-0.44484777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-0.96546121) q[1];
rz(-2.0270258) q[3];
sx q[3];
rz(-0.8884512) q[3];
sx q[3];
rz(1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(2.1202309) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9766831) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(3.0747736) q[0];
rz(1.5049997) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(0.79236275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6639858) q[1];
sx q[1];
rz(-2.7684282) q[1];
sx q[1];
rz(-0.54877703) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93123575) q[3];
sx q[3];
rz(-2.3124472) q[3];
sx q[3];
rz(-0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-0.41346082) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-2.5543509) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2205692) q[0];
sx q[0];
rz(-2.0221634) q[0];
sx q[0];
rz(1.736182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7325399) q[2];
sx q[2];
rz(-1.6502893) q[2];
sx q[2];
rz(0.0083991945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8752012) q[1];
sx q[1];
rz(-0.57386639) q[1];
sx q[1];
rz(0.89504524) q[1];
rz(-pi) q[2];
rz(2.5891853) q[3];
sx q[3];
rz(-0.93821412) q[3];
sx q[3];
rz(-1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(0.89909536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16222787) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(-2.2158932) q[0];
rz(-1.0663509) q[2];
sx q[2];
rz(-2.5214508) q[2];
sx q[2];
rz(-1.3401741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2832665) q[1];
sx q[1];
rz(-0.35052931) q[1];
sx q[1];
rz(1.9117029) q[1];
rz(-pi) q[2];
rz(-1.5702463) q[3];
sx q[3];
rz(-2.7714202) q[3];
sx q[3];
rz(-2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(-2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1966232) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(1.0265795) q[0];
x q[1];
rz(-2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(-3.0416833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(1.3089048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5554572) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(1.6554228) q[2];
sx q[2];
rz(-1.6783236) q[2];
sx q[2];
rz(2.8477737) q[2];
rz(0.11937033) q[3];
sx q[3];
rz(-1.3138249) q[3];
sx q[3];
rz(-1.8520595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];