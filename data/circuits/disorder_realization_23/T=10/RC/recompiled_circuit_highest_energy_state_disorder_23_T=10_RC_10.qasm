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
rz(0.5956369) q[0];
sx q[0];
rz(3.5516153) q[0];
sx q[0];
rz(8.6300996) q[0];
rz(2.088264) q[1];
sx q[1];
rz(-2.1894708) q[1];
sx q[1];
rz(1.9917537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2993752) q[0];
sx q[0];
rz(-1.2769852) q[0];
sx q[0];
rz(-0.037328193) q[0];
rz(-1.8914521) q[2];
sx q[2];
rz(-2.3806664) q[2];
sx q[2];
rz(0.93283949) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.036630298) q[1];
sx q[1];
rz(-1.399002) q[1];
sx q[1];
rz(-2.7378332) q[1];
x q[2];
rz(1.3142228) q[3];
sx q[3];
rz(-2.024641) q[3];
sx q[3];
rz(2.3875349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(2.4993994) q[2];
rz(2.0726974) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.4890891) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(-0.51134837) q[0];
rz(-1.6343575) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(1.7587597) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(14/(3*pi)) q[0];
sx q[0];
rz(-1.9535747) q[0];
sx q[0];
rz(2.9188927) q[0];
rz(-2.2425425) q[2];
sx q[2];
rz(-1.8457638) q[2];
sx q[2];
rz(3.0877047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1630711) q[1];
sx q[1];
rz(-2.6252665) q[1];
sx q[1];
rz(0.90746812) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17937998) q[3];
sx q[3];
rz(-0.53032833) q[3];
sx q[3];
rz(1.3190003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6766659) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(1.0627221) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17387834) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(3.0318731) q[0];
rz(2.8166855) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(-0.5330162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1905381) q[0];
sx q[0];
rz(-1.2742654) q[0];
sx q[0];
rz(-1.426479) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7092878) q[2];
sx q[2];
rz(-0.057261618) q[2];
sx q[2];
rz(-0.5106411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91777201) q[1];
sx q[1];
rz(-0.9448565) q[1];
sx q[1];
rz(-0.84898941) q[1];
rz(-pi) q[2];
rz(-2.332052) q[3];
sx q[3];
rz(-1.24461) q[3];
sx q[3];
rz(-2.1128138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32536062) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(1.4556966) q[2];
rz(3.0560737) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(-0.94676179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7607255) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(1.5707312) q[1];
sx q[1];
rz(-1.1474362) q[1];
sx q[1];
rz(3.0991203) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0642424) q[0];
sx q[0];
rz(-1.1288861) q[0];
sx q[0];
rz(1.7741706) q[0];
rz(-pi) q[1];
rz(-3.1257079) q[2];
sx q[2];
rz(-1.6611036) q[2];
sx q[2];
rz(1.7866194) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5407357) q[1];
sx q[1];
rz(-1.7607948) q[1];
sx q[1];
rz(0.79278391) q[1];
x q[2];
rz(-1.3298395) q[3];
sx q[3];
rz(-2.3577839) q[3];
sx q[3];
rz(2.2026521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(-0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1112082) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(0.17330387) q[0];
rz(-2.2068842) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-2.8876143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36080632) q[0];
sx q[0];
rz(-2.3179231) q[0];
sx q[0];
rz(-0.15374462) q[0];
rz(-pi) q[1];
rz(-2.5475909) q[2];
sx q[2];
rz(-1.2662167) q[2];
sx q[2];
rz(2.0128743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3601927) q[1];
sx q[1];
rz(-1.5915697) q[1];
sx q[1];
rz(2.4482577) q[1];
x q[2];
rz(3.0766904) q[3];
sx q[3];
rz(-1.9577454) q[3];
sx q[3];
rz(-2.153819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0678593) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(1.7013288) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.22236) q[3];
sx q[3];
rz(-3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9192231) q[0];
sx q[0];
rz(-1.8310522) q[0];
sx q[0];
rz(2.4134912) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(-2.5808835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9084204) q[0];
sx q[0];
rz(-1.5058555) q[0];
sx q[0];
rz(-3.0149595) q[0];
x q[1];
rz(1.4523023) q[2];
sx q[2];
rz(-1.0241707) q[2];
sx q[2];
rz(2.7989863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47885671) q[1];
sx q[1];
rz(-0.4534035) q[1];
sx q[1];
rz(-2.0824361) q[1];
rz(-pi) q[2];
rz(0.097499121) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(-1.1600375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59549436) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(2.0687436) q[2];
rz(0.32235518) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(-2.7505007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748465) q[0];
sx q[0];
rz(-1.5006737) q[0];
sx q[0];
rz(2.7333976) q[0];
rz(-2.8152668) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(-1.4761285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53083071) q[0];
sx q[0];
rz(-2.9617887) q[0];
sx q[0];
rz(1.7925949) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2017952) q[2];
sx q[2];
rz(-0.59518669) q[2];
sx q[2];
rz(-3.0734504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.020425107) q[1];
sx q[1];
rz(-1.5161991) q[1];
sx q[1];
rz(1.4538641) q[1];
rz(-pi) q[2];
rz(1.7404775) q[3];
sx q[3];
rz(-1.7998003) q[3];
sx q[3];
rz(1.218971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8866715) q[2];
sx q[2];
rz(-1.485606) q[2];
sx q[2];
rz(2.2754748) q[2];
rz(2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(-2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199468) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.2787) q[0];
rz(2.1389029) q[1];
sx q[1];
rz(-1.2607144) q[1];
sx q[1];
rz(-2.6967646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072413) q[0];
sx q[0];
rz(-1.0486516) q[0];
sx q[0];
rz(0.17052167) q[0];
rz(2.629822) q[2];
sx q[2];
rz(-2.3507042) q[2];
sx q[2];
rz(-2.1573586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5096566) q[1];
sx q[1];
rz(-2.3705814) q[1];
sx q[1];
rz(1.7421175) q[1];
rz(1.3108389) q[3];
sx q[3];
rz(-1.204986) q[3];
sx q[3];
rz(1.2838252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.727641) q[2];
sx q[2];
rz(0.50602305) q[2];
rz(2.1972726) q[3];
sx q[3];
rz(-0.025886141) q[3];
sx q[3];
rz(2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29189062) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(-3.131409) q[0];
rz(0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(0.082286509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5389494) q[0];
sx q[0];
rz(-1.2017439) q[0];
sx q[0];
rz(-2.6239388) q[0];
x q[1];
rz(1.2882223) q[2];
sx q[2];
rz(-2.4420441) q[2];
sx q[2];
rz(-2.4440773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7499867) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(2.8620666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.786834) q[3];
sx q[3];
rz(-0.91751614) q[3];
sx q[3];
rz(0.21524425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33132195) q[2];
sx q[2];
rz(-1.6931345) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(1.4823312) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3191147) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(2.5122232) q[0];
rz(1.7259701) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.9633044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4067766) q[0];
sx q[0];
rz(-1.1685851) q[0];
sx q[0];
rz(0.91535416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1599102) q[2];
sx q[2];
rz(-2.5576563) q[2];
sx q[2];
rz(0.46729014) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32543031) q[1];
sx q[1];
rz(-2.3419082) q[1];
sx q[1];
rz(-0.15284784) q[1];
rz(-pi) q[2];
rz(2.4201336) q[3];
sx q[3];
rz(-2.0155523) q[3];
sx q[3];
rz(-3.0196684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(-0.81864041) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(-2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4928987) q[0];
sx q[0];
rz(-2.0317827) q[0];
sx q[0];
rz(-1.6775525) q[0];
rz(1.875444) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(0.89546236) q[2];
sx q[2];
rz(-1.7747028) q[2];
sx q[2];
rz(1.5595421) q[2];
rz(2.9569382) q[3];
sx q[3];
rz(-1.481537) q[3];
sx q[3];
rz(2.8164345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
