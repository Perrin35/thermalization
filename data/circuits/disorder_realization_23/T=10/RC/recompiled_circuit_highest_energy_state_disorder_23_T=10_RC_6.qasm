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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(-2.3469143) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73939378) q[0];
sx q[0];
rz(-1.6065242) q[0];
sx q[0];
rz(1.2767919) q[0];
rz(-pi) q[1];
rz(0.29157421) q[2];
sx q[2];
rz(-0.85735029) q[2];
sx q[2];
rz(2.6387362) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.036630298) q[1];
sx q[1];
rz(-1.399002) q[1];
sx q[1];
rz(0.40375945) q[1];
rz(-pi) q[2];
rz(-0.46709316) q[3];
sx q[3];
rz(-1.8008999) q[3];
sx q[3];
rz(2.439374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6940234) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(0.64219323) q[2];
rz(1.0688952) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(-1.6525035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(-0.51134837) q[0];
rz(1.5072352) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(1.3828329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6561474) q[0];
sx q[0];
rz(-1.9535747) q[0];
sx q[0];
rz(0.22269999) q[0];
rz(-pi) q[1];
rz(-1.1452008) q[2];
sx q[2];
rz(-0.71768242) q[2];
sx q[2];
rz(1.8459783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1369776) q[1];
sx q[1];
rz(-1.8796693) q[1];
sx q[1];
rz(-1.1501794) q[1];
x q[2];
rz(2.6182943) q[3];
sx q[3];
rz(-1.4804258) q[3];
sx q[3];
rz(-0.096631526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4649268) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(1.0627221) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(-3.0265871) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17387834) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(-3.0318731) q[0];
rz(2.8166855) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(0.5330162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.479411) q[0];
sx q[0];
rz(-1.7087738) q[0];
sx q[0];
rz(0.29946391) q[0];
rz(-pi) q[1];
rz(1.5948086) q[2];
sx q[2];
rz(-1.622785) q[2];
sx q[2];
rz(2.1980224) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.24025) q[1];
sx q[1];
rz(-2.2247215) q[1];
sx q[1];
rz(-2.4010977) q[1];
rz(-pi) q[2];
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
rz(2.816232) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(-3.0560737) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(0.94676179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
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
rz(2.6149276) q[0];
sx q[0];
rz(-0.4836429) q[0];
sx q[0];
rz(2.7381104) q[0];
rz(-pi) q[1];
rz(-1.7444403) q[2];
sx q[2];
rz(-0.091689907) q[2];
sx q[2];
rz(1.1806115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9876696) q[1];
sx q[1];
rz(-2.3312285) q[1];
sx q[1];
rz(-0.26370476) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3399502) q[3];
sx q[3];
rz(-1.7400734) q[3];
sx q[3];
rz(-2.3374286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2461569) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(2.4401276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0303845) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(-2.9682888) q[0];
rz(-2.2068842) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-2.8876143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58502176) q[0];
sx q[0];
rz(-0.75977548) q[0];
sx q[0];
rz(-1.7346481) q[0];
rz(1.9333657) q[2];
sx q[2];
rz(-2.1340279) q[2];
sx q[2];
rz(-2.8993487) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9137302) q[1];
sx q[1];
rz(-2.263952) q[1];
sx q[1];
rz(1.597803) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0766904) q[3];
sx q[3];
rz(-1.1838473) q[3];
sx q[3];
rz(0.98777366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(0.67617792) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.22236958) q[0];
sx q[0];
rz(-1.8310522) q[0];
sx q[0];
rz(-0.72810143) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-0.89465061) q[1];
sx q[1];
rz(2.5808835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9084204) q[0];
sx q[0];
rz(-1.5058555) q[0];
sx q[0];
rz(3.0149595) q[0];
rz(-pi) q[1];
rz(-0.19189437) q[2];
sx q[2];
rz(-0.55804306) q[2];
sx q[2];
rz(0.11745889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5594031) q[1];
sx q[1];
rz(-1.786937) q[1];
sx q[1];
rz(1.1690421) q[1];
x q[2];
rz(0.097499121) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(1.9815552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59549436) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(2.0687436) q[2];
rz(-2.8192375) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(0.391092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748465) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(0.32632581) q[1];
sx q[1];
rz(-0.34714454) q[1];
sx q[1];
rz(-1.6654642) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583019) q[0];
sx q[0];
rz(-1.6101478) q[0];
sx q[0];
rz(1.7462867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0075241) q[2];
sx q[2];
rz(-1.7744229) q[2];
sx q[2];
rz(1.3290392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1154815) q[1];
sx q[1];
rz(-3.012595) q[1];
sx q[1];
rz(1.1327101) q[1];
x q[2];
rz(-2.514634) q[3];
sx q[3];
rz(-2.857465) q[3];
sx q[3];
rz(-1.8654902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8866715) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(-2.2754748) q[2];
rz(-2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199468) q[0];
sx q[0];
rz(-0.045504657) q[0];
sx q[0];
rz(1.8628927) q[0];
rz(-2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-2.6967646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507752) q[0];
sx q[0];
rz(-1.4231761) q[0];
sx q[0];
rz(-1.0423129) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1110531) q[2];
sx q[2];
rz(-0.90219775) q[2];
sx q[2];
rz(2.8313314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39522102) q[1];
sx q[1];
rz(-0.81391293) q[1];
sx q[1];
rz(-2.9774352) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8307537) q[3];
sx q[3];
rz(-1.204986) q[3];
sx q[3];
rz(1.8577675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7030846) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-2.6355696) q[2];
rz(-2.1972726) q[3];
sx q[3];
rz(-0.025886141) q[3];
sx q[3];
rz(-2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(-3.131409) q[0];
rz(0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-3.0593061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40348379) q[0];
sx q[0];
rz(-2.5157564) q[0];
sx q[0];
rz(2.4781669) q[0];
rz(-0.89110969) q[2];
sx q[2];
rz(-1.3902877) q[2];
sx q[2];
rz(-1.0918822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.391606) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(0.27952607) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27307053) q[3];
sx q[3];
rz(-0.6830754) q[3];
sx q[3];
rz(-0.13126743) q[3];
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
rz(-1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.822478) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(0.6293695) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(1.9633044) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4067766) q[0];
sx q[0];
rz(-1.1685851) q[0];
sx q[0];
rz(-2.2262385) q[0];
x q[1];
rz(2.0731931) q[2];
sx q[2];
rz(-1.2594688) q[2];
sx q[2];
rz(-2.5466998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7892742) q[1];
sx q[1];
rz(-1.4613918) q[1];
sx q[1];
rz(-2.3477702) q[1];
x q[2];
rz(-2.4201336) q[3];
sx q[3];
rz(-1.1260403) q[3];
sx q[3];
rz(0.12192425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(-2.3229522) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64869399) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(-1.875444) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(-0.89546236) q[2];
sx q[2];
rz(-1.3668899) q[2];
sx q[2];
rz(-1.5820506) q[2];
rz(-0.18465445) q[3];
sx q[3];
rz(-1.481537) q[3];
sx q[3];
rz(2.8164345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
