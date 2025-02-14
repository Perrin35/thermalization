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
rz(5.6018957) q[1];
sx q[1];
rz(8.3991474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5012623) q[0];
sx q[0];
rz(-0.20168951) q[0];
sx q[0];
rz(1.8224688) q[0];
rz(-pi) q[1];
rz(-2.9812987) q[2];
sx q[2];
rz(-2.829814) q[2];
sx q[2];
rz(2.3202053) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.611022) q[1];
sx q[1];
rz(-2.321302) q[1];
sx q[1];
rz(1.6802701) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4769931) q[3];
sx q[3];
rz(-2.0076111) q[3];
sx q[3];
rz(0.38540456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2166298) q[2];
sx q[2];
rz(-1.3684042) q[2];
sx q[2];
rz(1.6538357) q[2];
rz(-1.1612859) q[3];
sx q[3];
rz(-2.1372644) q[3];
sx q[3];
rz(2.0465093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.0733114) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(-2.708129) q[0];
rz(1.06217) q[1];
sx q[1];
rz(-1.5824317) q[1];
sx q[1];
rz(0.72057048) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67340785) q[0];
sx q[0];
rz(-1.6077198) q[0];
sx q[0];
rz(1.4440749) q[0];
rz(-2.8938204) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(1.2733851) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99834187) q[1];
sx q[1];
rz(-1.9228808) q[1];
sx q[1];
rz(-1.4125173) q[1];
rz(-pi) q[2];
rz(1.0834915) q[3];
sx q[3];
rz(-1.8326899) q[3];
sx q[3];
rz(2.0341121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34708193) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(-2.1573055) q[2];
rz(0.85683626) q[3];
sx q[3];
rz(-2.555116) q[3];
sx q[3];
rz(-1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(3.1318829) q[1];
sx q[1];
rz(-1.566889) q[1];
sx q[1];
rz(-0.28365338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6742636) q[0];
sx q[0];
rz(-1.0630075) q[0];
sx q[0];
rz(1.171815) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7314807) q[2];
sx q[2];
rz(-1.3734439) q[2];
sx q[2];
rz(1.6244217) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4176826) q[1];
sx q[1];
rz(-1.7387448) q[1];
sx q[1];
rz(-1.8803094) q[1];
x q[2];
rz(-1.4693442) q[3];
sx q[3];
rz(-1.5939004) q[3];
sx q[3];
rz(0.35233179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126605) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(1.5869045) q[0];
rz(-2.9117865) q[1];
sx q[1];
rz(-0.72139144) q[1];
sx q[1];
rz(1.302964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763846) q[0];
sx q[0];
rz(-2.1094649) q[0];
sx q[0];
rz(-1.046565) q[0];
x q[1];
rz(-1.8660061) q[2];
sx q[2];
rz(-1.3129932) q[2];
sx q[2];
rz(0.38391963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.069217056) q[1];
sx q[1];
rz(-1.0472968) q[1];
sx q[1];
rz(1.545946) q[1];
x q[2];
rz(-2.0992364) q[3];
sx q[3];
rz(-2.5109249) q[3];
sx q[3];
rz(-1.3321277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0064696781) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(2.6378677) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(-1.5776207) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85394323) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(-0.26643878) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(2.0727167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7212413) q[0];
sx q[0];
rz(-1.9801894) q[0];
sx q[0];
rz(-0.73914011) q[0];
x q[1];
rz(-1.8891625) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(-0.81046644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68764781) q[1];
sx q[1];
rz(-1.3140716) q[1];
sx q[1];
rz(-2.686108) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74917082) q[3];
sx q[3];
rz(-0.35532829) q[3];
sx q[3];
rz(-2.4859249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2499007) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(-1.8788667) q[3];
sx q[3];
rz(-1.6014674) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1332909) q[0];
sx q[0];
rz(-2.1228078) q[0];
sx q[0];
rz(-0.89023501) q[0];
rz(-1.4324191) q[1];
sx q[1];
rz(-1.018367) q[1];
sx q[1];
rz(-0.11996809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16968834) q[0];
sx q[0];
rz(-1.7154046) q[0];
sx q[0];
rz(0.55938104) q[0];
x q[1];
rz(-0.34949686) q[2];
sx q[2];
rz(-1.9781402) q[2];
sx q[2];
rz(-1.1098338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90112061) q[1];
sx q[1];
rz(-1.3252621) q[1];
sx q[1];
rz(-0.66048173) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7637353) q[3];
sx q[3];
rz(-1.0385374) q[3];
sx q[3];
rz(1.618945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.031781901) q[2];
sx q[2];
rz(-1.4795156) q[2];
sx q[2];
rz(-0.39153448) q[2];
rz(1.7560962) q[3];
sx q[3];
rz(-0.23844312) q[3];
sx q[3];
rz(-2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(1.3054003) q[0];
rz(-2.6370559) q[1];
sx q[1];
rz(-1.4089818) q[1];
sx q[1];
rz(0.22444589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2674057) q[0];
sx q[0];
rz(-1.5541557) q[0];
sx q[0];
rz(-1.5550625) q[0];
x q[1];
rz(-1.7845909) q[2];
sx q[2];
rz(-2.0540621) q[2];
sx q[2];
rz(-2.4501806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45221165) q[1];
sx q[1];
rz(-1.2465917) q[1];
sx q[1];
rz(0.28181847) q[1];
rz(1.2013632) q[3];
sx q[3];
rz(-1.2681586) q[3];
sx q[3];
rz(1.3898046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-0.9474729) q[2];
rz(0.13253658) q[3];
sx q[3];
rz(-0.65968502) q[3];
sx q[3];
rz(-0.98602492) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9082311) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(2.76873) q[0];
rz(-2.8064959) q[1];
sx q[1];
rz(-1.2058039) q[1];
sx q[1];
rz(-2.0705409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.459254) q[0];
sx q[0];
rz(-1.5594578) q[0];
sx q[0];
rz(-3.0724694) q[0];
x q[1];
rz(-0.33800563) q[2];
sx q[2];
rz(-0.57799423) q[2];
sx q[2];
rz(-1.4678616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4651686) q[1];
sx q[1];
rz(-1.3260726) q[1];
sx q[1];
rz(1.8108435) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0464977) q[3];
sx q[3];
rz(-1.250306) q[3];
sx q[3];
rz(1.2456428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7143453) q[2];
sx q[2];
rz(-1.599396) q[2];
sx q[2];
rz(-1.0672807) q[2];
rz(0.37211564) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(0.20017643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9099971) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.7344612) q[0];
rz(1.0435957) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(-2.6503906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022302786) q[0];
sx q[0];
rz(-0.91141846) q[0];
sx q[0];
rz(-2.7727866) q[0];
x q[1];
rz(0.78989002) q[2];
sx q[2];
rz(-2.9040376) q[2];
sx q[2];
rz(-1.6611163) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8206827) q[1];
sx q[1];
rz(-1.2315799) q[1];
sx q[1];
rz(1.7103819) q[1];
x q[2];
rz(-1.7513428) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(-2.2425687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4444943) q[2];
sx q[2];
rz(-1.631087) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(0.77110428) q[3];
sx q[3];
rz(-1.6090569) q[3];
sx q[3];
rz(0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20030178) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(-1.5197165) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(2.4470952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8246758) q[0];
sx q[0];
rz(-1.7986713) q[0];
sx q[0];
rz(0.014046305) q[0];
x q[1];
rz(2.6171309) q[2];
sx q[2];
rz(-1.6860644) q[2];
sx q[2];
rz(2.070321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81238231) q[1];
sx q[1];
rz(-0.33420104) q[1];
sx q[1];
rz(-1.5370395) q[1];
x q[2];
rz(-1.9574088) q[3];
sx q[3];
rz(-2.2848633) q[3];
sx q[3];
rz(0.82502818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0684315) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(0.51363242) q[2];
rz(0.3178151) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(-1.0191127) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014872) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(1.2263251) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(0.60811715) q[2];
sx q[2];
rz(-0.79610745) q[2];
sx q[2];
rz(2.983762) q[2];
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
