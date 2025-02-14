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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(0.96487784) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2478188) q[0];
sx q[0];
rz(-2.4501094) q[0];
sx q[0];
rz(-1.0214424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9994316) q[2];
sx q[2];
rz(-1.8920676) q[2];
sx q[2];
rz(0.69883332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4360848) q[1];
sx q[1];
rz(-0.75172808) q[1];
sx q[1];
rz(-0.99940325) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0003188) q[3];
sx q[3];
rz(-2.5577099) q[3];
sx q[3];
rz(-0.48030765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45399484) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(0.54443693) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(2.4936567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(-2.6065705) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-0.63553634) q[1];
sx q[1];
rz(1.7417057) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28093058) q[0];
sx q[0];
rz(-1.1854404) q[0];
sx q[0];
rz(2.0549261) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1715047) q[2];
sx q[2];
rz(-1.0905438) q[2];
sx q[2];
rz(2.9532972) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4606159) q[1];
sx q[1];
rz(-1.9505122) q[1];
sx q[1];
rz(1.3028076) q[1];
rz(1.0008282) q[3];
sx q[3];
rz(-2.1080351) q[3];
sx q[3];
rz(-1.9443823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(-0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(-1.6287623) q[0];
rz(-0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(-2.0294752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84217269) q[0];
sx q[0];
rz(-1.5779691) q[0];
sx q[0];
rz(0.0070863574) q[0];
rz(-2.1826964) q[2];
sx q[2];
rz(-3.009353) q[2];
sx q[2];
rz(2.7041777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97858799) q[1];
sx q[1];
rz(-1.7573865) q[1];
sx q[1];
rz(-0.48689894) q[1];
rz(2.882631) q[3];
sx q[3];
rz(-1.9230712) q[3];
sx q[3];
rz(1.5332298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(-0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(-0.68921047) q[0];
rz(-2.8269732) q[1];
sx q[1];
rz(-2.2210821) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5310881) q[0];
sx q[0];
rz(-2.304739) q[0];
sx q[0];
rz(0.42829163) q[0];
x q[1];
rz(-2.7233549) q[2];
sx q[2];
rz(-2.9091638) q[2];
sx q[2];
rz(1.7247891) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8103579) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-0.15366252) q[1];
rz(-1.9288428) q[3];
sx q[3];
rz(-2.6302377) q[3];
sx q[3];
rz(2.8062888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(-2.5856384) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81412643) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(0.59447527) q[0];
rz(2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-2.1655703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291479) q[0];
sx q[0];
rz(-1.2093822) q[0];
sx q[0];
rz(2.0301719) q[0];
x q[1];
rz(-3.1107509) q[2];
sx q[2];
rz(-2.387306) q[2];
sx q[2];
rz(-0.13067836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4415293) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(1.1928012) q[1];
x q[2];
rz(-2.8401883) q[3];
sx q[3];
rz(-1.7782974) q[3];
sx q[3];
rz(1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(2.5308934) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.82764757) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(1.4661283) q[0];
rz(0.77955359) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-0.73878845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2125428) q[0];
sx q[0];
rz(-1.7746141) q[0];
sx q[0];
rz(-1.8575559) q[0];
rz(-pi) q[1];
x q[1];
rz(0.03348695) q[2];
sx q[2];
rz(-2.8346363) q[2];
sx q[2];
rz(0.54881664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9924817) q[1];
sx q[1];
rz(-2.457239) q[1];
sx q[1];
rz(-2.1696287) q[1];
rz(-pi) q[2];
rz(0.8035369) q[3];
sx q[3];
rz(-2.6430686) q[3];
sx q[3];
rz(-0.62832181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61591992) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(2.681813) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.1427243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2343242) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(-2.1211076) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(0.78757706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575229) q[0];
sx q[0];
rz(-2.2057187) q[0];
sx q[0];
rz(1.8674839) q[0];
rz(-pi) q[1];
rz(1.4992981) q[2];
sx q[2];
rz(-1.8282969) q[2];
sx q[2];
rz(-1.6704287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5912093) q[1];
sx q[1];
rz(-1.0799284) q[1];
sx q[1];
rz(0.20259133) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1586886) q[3];
sx q[3];
rz(-1.5421151) q[3];
sx q[3];
rz(0.3405638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3202177) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(2.5569432) q[2];
rz(2.4892877) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(1.339213) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(-0.32824326) q[0];
rz(2.7499061) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(-1.4788871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680964) q[0];
sx q[0];
rz(-2.7910821) q[0];
sx q[0];
rz(-2.4309733) q[0];
rz(-2.4166862) q[2];
sx q[2];
rz(-2.4838243) q[2];
sx q[2];
rz(0.36054128) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6819671) q[1];
sx q[1];
rz(-0.29469583) q[1];
sx q[1];
rz(1.2380283) q[1];
rz(-1.9401266) q[3];
sx q[3];
rz(-1.7672667) q[3];
sx q[3];
rz(-0.35614355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0216003) q[2];
sx q[2];
rz(-1.768521) q[2];
sx q[2];
rz(0.78869406) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(0.23421639) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(0.73807565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2797425) q[0];
sx q[0];
rz(-1.8075917) q[0];
sx q[0];
rz(-2.7351232) q[0];
rz(-pi) q[1];
rz(2.5665087) q[2];
sx q[2];
rz(-0.7204537) q[2];
sx q[2];
rz(1.6937814) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.187584) q[1];
sx q[1];
rz(-0.52032797) q[1];
sx q[1];
rz(-2.1854758) q[1];
x q[2];
rz(-1.290019) q[3];
sx q[3];
rz(-2.2178136) q[3];
sx q[3];
rz(0.10296497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(1.8812995) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(-0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(2.2743478) q[0];
rz(-2.1278837) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(-1.0702466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1742221) q[0];
sx q[0];
rz(-1.692932) q[0];
sx q[0];
rz(-1.4184667) q[0];
x q[1];
rz(-2.321601) q[2];
sx q[2];
rz(-2.242616) q[2];
sx q[2];
rz(-1.9209678) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5336736) q[1];
sx q[1];
rz(-1.1864788) q[1];
sx q[1];
rz(2.584105) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66859122) q[3];
sx q[3];
rz(-0.1771268) q[3];
sx q[3];
rz(-1.9884381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16537198) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(1.9099859) q[2];
rz(1.6407137) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-2.7267743) q[1];
sx q[1];
rz(-1.7638313) q[1];
sx q[1];
rz(1.5324963) q[1];
rz(2.6957569) q[2];
sx q[2];
rz(-1.5210963) q[2];
sx q[2];
rz(-1.7901909) q[2];
rz(-2.5255193) q[3];
sx q[3];
rz(-0.85312927) q[3];
sx q[3];
rz(2.7240172) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
