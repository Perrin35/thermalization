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
rz(-2.73157) q[0];
sx q[0];
rz(-0.79467839) q[0];
rz(2.088264) q[1];
sx q[1];
rz(-2.1894708) q[1];
sx q[1];
rz(-1.149839) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4276221) q[0];
sx q[0];
rz(-0.29610482) q[0];
sx q[0];
rz(1.6935191) q[0];
rz(-pi) q[1];
rz(2.8500184) q[2];
sx q[2];
rz(-0.85735029) q[2];
sx q[2];
rz(-2.6387362) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1049624) q[1];
sx q[1];
rz(-1.399002) q[1];
sx q[1];
rz(0.40375945) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46709316) q[3];
sx q[3];
rz(-1.3406928) q[3];
sx q[3];
rz(-2.439374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(-2.4993994) q[2];
rz(2.0726974) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(-1.4890891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25817961) q[0];
sx q[0];
rz(-3.1386107) q[0];
sx q[0];
rz(-2.6302443) q[0];
rz(1.6343575) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(1.7587597) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6561474) q[0];
sx q[0];
rz(-1.1880179) q[0];
sx q[0];
rz(2.9188927) q[0];
x q[1];
rz(2.7956656) q[2];
sx q[2];
rz(-2.2130161) q[2];
sx q[2];
rz(-1.3042892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1369776) q[1];
sx q[1];
rz(-1.2619234) q[1];
sx q[1];
rz(1.9914133) q[1];
x q[2];
rz(2.9622127) q[3];
sx q[3];
rz(-2.6112643) q[3];
sx q[3];
rz(-1.3190003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6766659) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(3.0265871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677143) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(-3.0318731) q[0];
rz(2.8166855) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(-0.5330162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66218161) q[0];
sx q[0];
rz(-1.4328188) q[0];
sx q[0];
rz(-2.8421287) q[0];
x q[1];
rz(-1.5948086) q[2];
sx q[2];
rz(-1.5188076) q[2];
sx q[2];
rz(-0.94357027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9013427) q[1];
sx q[1];
rz(-0.91687119) q[1];
sx q[1];
rz(-0.74049495) q[1];
x q[2];
rz(0.80954062) q[3];
sx q[3];
rz(-1.8969826) q[3];
sx q[3];
rz(-1.0287788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32536062) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(-3.0560737) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(2.1948309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38086712) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(-0.72232676) q[0];
rz(-1.5708615) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-3.0991203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149276) q[0];
sx q[0];
rz(-2.6579497) q[0];
sx q[0];
rz(2.7381104) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4804777) q[2];
sx q[2];
rz(-1.5866163) q[2];
sx q[2];
rz(2.9243369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60085697) q[1];
sx q[1];
rz(-1.3807978) q[1];
sx q[1];
rz(-2.3488087) q[1];
rz(0.23353429) q[3];
sx q[3];
rz(-0.81538768) q[3];
sx q[3];
rz(2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8954358) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(2.4401276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0303845) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(-0.17330387) q[0];
rz(2.2068842) q[1];
sx q[1];
rz(-0.84988958) q[1];
sx q[1];
rz(-2.8876143) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7807863) q[0];
sx q[0];
rz(-0.82366952) q[0];
sx q[0];
rz(-0.15374462) q[0];
rz(-pi) q[1];
rz(2.5475909) q[2];
sx q[2];
rz(-1.8753759) q[2];
sx q[2];
rz(2.0128743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9137302) q[1];
sx q[1];
rz(-0.87764064) q[1];
sx q[1];
rz(-1.597803) q[1];
rz(0.06490223) q[3];
sx q[3];
rz(-1.1838473) q[3];
sx q[3];
rz(0.98777366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(1.7013288) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(-0.079782709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22236958) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(-2.4134912) q[0];
rz(-1.9758196) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(2.5808835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8122305) q[0];
sx q[0];
rz(-1.6971611) q[0];
sx q[0];
rz(1.6362599) q[0];
rz(-1.4523023) q[2];
sx q[2];
rz(-2.117422) q[2];
sx q[2];
rz(2.7989863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47885671) q[1];
sx q[1];
rz(-2.6881892) q[1];
sx q[1];
rz(-1.0591566) q[1];
rz(-pi) q[2];
rz(-3.0440935) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(-1.1600375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.76674616) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(0.40819502) q[0];
rz(2.8152668) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(-1.6654642) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583019) q[0];
sx q[0];
rz(-1.5314449) q[0];
sx q[0];
rz(-1.7462867) q[0];
rz(1.0075241) q[2];
sx q[2];
rz(-1.3671698) q[2];
sx q[2];
rz(-1.3290392) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0261111) q[1];
sx q[1];
rz(-3.012595) q[1];
sx q[1];
rz(-2.0088825) q[1];
rz(1.4011151) q[3];
sx q[3];
rz(-1.3417923) q[3];
sx q[3];
rz(-1.9226216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8866715) q[2];
sx q[2];
rz(-1.6559867) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12164584) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.8628927) q[0];
rz(-2.1389029) q[1];
sx q[1];
rz(-1.2607144) q[1];
sx q[1];
rz(2.6967646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76679517) q[0];
sx q[0];
rz(-0.54682362) q[0];
sx q[0];
rz(-1.8575791) q[0];
x q[1];
rz(1.1110531) q[2];
sx q[2];
rz(-2.2393949) q[2];
sx q[2];
rz(-0.3102613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0792744) q[1];
sx q[1];
rz(-1.6898815) q[1];
sx q[1];
rz(-2.3344385) q[1];
rz(0.59101848) q[3];
sx q[3];
rz(-0.44535397) q[3];
sx q[3];
rz(1.9231918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43850809) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(0.50602305) q[2];
rz(2.1972726) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29189062) q[0];
sx q[0];
rz(-2.1493752) q[0];
sx q[0];
rz(0.010183656) q[0];
rz(0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-3.0593061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40348379) q[0];
sx q[0];
rz(-0.62583621) q[0];
sx q[0];
rz(-2.4781669) q[0];
rz(-0.89110969) q[2];
sx q[2];
rz(-1.3902877) q[2];
sx q[2];
rz(-1.0918822) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6895143) q[1];
sx q[1];
rz(-1.6321215) q[1];
sx q[1];
rz(0.21558) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8685221) q[3];
sx q[3];
rz(-0.6830754) q[3];
sx q[3];
rz(-3.0103252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8102707) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(-1.4823312) q[3];
sx q[3];
rz(-0.63153657) q[3];
sx q[3];
rz(0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3191147) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(1.9633044) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4067766) q[0];
sx q[0];
rz(-1.1685851) q[0];
sx q[0];
rz(0.91535416) q[0];
rz(-0.35188108) q[2];
sx q[2];
rz(-1.0946254) q[2];
sx q[2];
rz(1.9989524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32543031) q[1];
sx q[1];
rz(-0.79968444) q[1];
sx q[1];
rz(-0.15284784) q[1];
rz(-2.516516) q[3];
sx q[3];
rz(-2.3155594) q[3];
sx q[3];
rz(0.99398289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(-2.3229522) q[2];
rz(1.5107752) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(-0.40217933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(-1.2661487) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(-1.2513592) q[2];
sx q[2];
rz(-0.70079679) q[2];
sx q[2];
rz(0.23637017) q[2];
rz(-2.6880445) q[3];
sx q[3];
rz(-2.9367179) q[3];
sx q[3];
rz(-2.3412326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
