OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(-2.4617221) q[0];
sx q[0];
rz(0.1006861) q[0];
rz(0.45473948) q[1];
sx q[1];
rz(-0.68128959) q[1];
sx q[1];
rz(-1.0256306) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89698638) q[0];
sx q[0];
rz(-1.375544) q[0];
sx q[0];
rz(3.0907187) q[0];
rz(-pi) q[1];
rz(1.5194015) q[2];
sx q[2];
rz(-1.2631466) q[2];
sx q[2];
rz(-2.4884698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.034596171) q[1];
sx q[1];
rz(-1.6507848) q[1];
sx q[1];
rz(-2.3880915) q[1];
x q[2];
rz(2.493606) q[3];
sx q[3];
rz(-0.77666908) q[3];
sx q[3];
rz(-2.4512993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92496282) q[2];
sx q[2];
rz(-1.3684042) q[2];
sx q[2];
rz(1.6538357) q[2];
rz(-1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(2.0465093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733114) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(2.708129) q[0];
rz(1.06217) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(-0.72057048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61536371) q[0];
sx q[0];
rz(-3.0096292) q[0];
sx q[0];
rz(-1.855164) q[0];
x q[1];
rz(2.1042688) q[2];
sx q[2];
rz(-1.7852655) q[2];
sx q[2];
rz(2.7191424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99834187) q[1];
sx q[1];
rz(-1.2187119) q[1];
sx q[1];
rz(1.4125173) q[1];
x q[2];
rz(-2.8470542) q[3];
sx q[3];
rz(-1.1014767) q[3];
sx q[3];
rz(-0.32696262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7945107) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(-0.98428717) q[2];
rz(-0.85683626) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(1.8872567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03265753) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6742636) q[0];
sx q[0];
rz(-2.0785851) q[0];
sx q[0];
rz(-1.9697777) q[0];
rz(-0.1998603) q[2];
sx q[2];
rz(-1.7283354) q[2];
sx q[2];
rz(-0.021856088) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3481118) q[1];
sx q[1];
rz(-1.2657796) q[1];
sx q[1];
rz(-2.9654337) q[1];
rz(-pi) q[2];
rz(1.795122) q[3];
sx q[3];
rz(-3.0375518) q[3];
sx q[3];
rz(-1.4416172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11188406) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(-0.65579826) q[2];
rz(2.9295975) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(-1.9688212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(1.5869045) q[0];
rz(0.2298062) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.302964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36520805) q[0];
sx q[0];
rz(-1.0321277) q[0];
sx q[0];
rz(2.0950277) q[0];
rz(-pi) q[1];
rz(-1.2755865) q[2];
sx q[2];
rz(-1.8285995) q[2];
sx q[2];
rz(0.38391963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.069217056) q[1];
sx q[1];
rz(-1.0472968) q[1];
sx q[1];
rz(-1.5956466) q[1];
rz(-pi) q[2];
rz(-1.0082208) q[3];
sx q[3];
rz(-1.2689212) q[3];
sx q[3];
rz(2.9396597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.135123) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(-2.6378677) q[2];
rz(-1.5375562) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(-1.5776207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.85394323) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(2.8751539) q[1];
sx q[1];
rz(-0.8100422) q[1];
sx q[1];
rz(-2.0727167) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3392068) q[0];
sx q[0];
rz(-0.9046208) q[0];
sx q[0];
rz(-1.0399014) q[0];
x q[1];
rz(-1.2524302) q[2];
sx q[2];
rz(-2.3471544) q[2];
sx q[2];
rz(2.3311262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4539448) q[1];
sx q[1];
rz(-1.8275211) q[1];
sx q[1];
rz(0.4554847) q[1];
rz(-pi) q[2];
rz(-2.3924218) q[3];
sx q[3];
rz(-2.7862644) q[3];
sx q[3];
rz(-0.65566777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8916919) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(2.967584) q[2];
rz(1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(1.5866097) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1332909) q[0];
sx q[0];
rz(-2.1228078) q[0];
sx q[0];
rz(-2.2513576) q[0];
rz(1.7091735) q[1];
sx q[1];
rz(-1.018367) q[1];
sx q[1];
rz(-0.11996809) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1748808) q[0];
sx q[0];
rz(-2.565756) q[0];
sx q[0];
rz(2.8737646) q[0];
x q[1];
rz(-1.1402848) q[2];
sx q[2];
rz(-1.8906279) q[2];
sx q[2];
rz(-2.8240273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6585582) q[1];
sx q[1];
rz(-0.93340988) q[1];
sx q[1];
rz(-1.263522) q[1];
rz(-pi) q[2];
rz(0.54048583) q[3];
sx q[3];
rz(-1.4048164) q[3];
sx q[3];
rz(0.14697187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1098108) q[2];
sx q[2];
rz(-1.4795156) q[2];
sx q[2];
rz(-2.7500582) q[2];
rz(1.7560962) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(-0.63867205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93381298) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.8361924) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(2.9171468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674057) q[0];
sx q[0];
rz(-1.5541557) q[0];
sx q[0];
rz(1.5550625) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6487892) q[2];
sx q[2];
rz(-1.3818008) q[2];
sx q[2];
rz(0.97993055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9516884) q[1];
sx q[1];
rz(-2.7153141) q[1];
sx q[1];
rz(2.262102) q[1];
x q[2];
rz(-1.9402294) q[3];
sx q[3];
rz(-1.2681586) q[3];
sx q[3];
rz(1.3898046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-0.9474729) q[2];
rz(-3.0090561) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(-2.1555677) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2333616) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(0.37286266) q[0];
rz(-0.33509675) q[1];
sx q[1];
rz(-1.2058039) q[1];
sx q[1];
rz(2.0705409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27386943) q[0];
sx q[0];
rz(-0.070045538) q[0];
sx q[0];
rz(-2.9788736) q[0];
rz(0.55166371) q[2];
sx q[2];
rz(-1.3886189) q[2];
sx q[2];
rz(-0.38924402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4651686) q[1];
sx q[1];
rz(-1.3260726) q[1];
sx q[1];
rz(1.3307491) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8926431) q[3];
sx q[3];
rz(-1.4805571) q[3];
sx q[3];
rz(-2.7863996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7143453) q[2];
sx q[2];
rz(-1.599396) q[2];
sx q[2];
rz(-1.0672807) q[2];
rz(2.769477) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(-0.20017643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2315955) q[0];
sx q[0];
rz(-2.5439926) q[0];
sx q[0];
rz(1.4071314) q[0];
rz(-2.097997) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(2.6503906) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564767) q[0];
sx q[0];
rz(-2.399754) q[0];
sx q[0];
rz(-2.0061352) q[0];
x q[1];
rz(-0.16881659) q[2];
sx q[2];
rz(-1.7387317) q[2];
sx q[2];
rz(2.4562648) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.421153) q[1];
sx q[1];
rz(-2.7758088) q[1];
sx q[1];
rz(0.37558742) q[1];
rz(-pi) q[2];
rz(1.7513428) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(-0.89902395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6970984) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(-0.77110428) q[3];
sx q[3];
rz(-1.6090569) q[3];
sx q[3];
rz(-0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20030178) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(-1.5197165) q[0];
rz(-1.8408076) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(-0.69449743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8845399) q[0];
sx q[0];
rz(-1.5844795) q[0];
sx q[0];
rz(1.3428997) q[0];
x q[1];
rz(-1.4378261) q[2];
sx q[2];
rz(-2.0914234) q[2];
sx q[2];
rz(0.43308738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.84811454) q[1];
sx q[1];
rz(-1.2367931) q[1];
sx q[1];
rz(3.1298742) q[1];
rz(-pi) q[2];
rz(1.1841838) q[3];
sx q[3];
rz(-2.2848633) q[3];
sx q[3];
rz(0.82502818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0731611) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(-0.51363242) q[2];
rz(-2.8237776) q[3];
sx q[3];
rz(-1.5376667) q[3];
sx q[3];
rz(-2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.9152676) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(2.4438159) q[2];
sx q[2];
rz(-1.991376) q[2];
sx q[2];
rz(-2.1817653) q[2];
rz(-1.574434) q[3];
sx q[3];
rz(-2.6156577) q[3];
sx q[3];
rz(-1.2729264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
