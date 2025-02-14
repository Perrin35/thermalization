OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98193869) q[0];
sx q[0];
rz(-1.5768134) q[0];
sx q[0];
rz(-1.2793581) q[0];
rz(0.98969412) q[1];
sx q[1];
rz(-1.1969748) q[1];
sx q[1];
rz(2.9762414) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26583126) q[0];
sx q[0];
rz(-0.94801694) q[0];
sx q[0];
rz(2.7490007) q[0];
rz(-2.697359) q[2];
sx q[2];
rz(-1.1886676) q[2];
sx q[2];
rz(-1.6319147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0667116) q[1];
sx q[1];
rz(-2.0600658) q[1];
sx q[1];
rz(2.679945) q[1];
rz(1.8824485) q[3];
sx q[3];
rz(-1.3247648) q[3];
sx q[3];
rz(-0.450799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87236658) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(-1.9744138) q[2];
rz(0.40927467) q[3];
sx q[3];
rz(-2.8661178) q[3];
sx q[3];
rz(-0.95612139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072227612) q[0];
sx q[0];
rz(-0.15412155) q[0];
sx q[0];
rz(1.7896205) q[0];
rz(1.6701472) q[1];
sx q[1];
rz(-0.92266005) q[1];
sx q[1];
rz(-0.89554375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2038747) q[0];
sx q[0];
rz(-0.26709891) q[0];
sx q[0];
rz(3.1338723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15094916) q[2];
sx q[2];
rz(-0.47892919) q[2];
sx q[2];
rz(-0.72266173) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18691639) q[1];
sx q[1];
rz(-2.3136931) q[1];
sx q[1];
rz(2.1174341) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9905521) q[3];
sx q[3];
rz(-1.7872056) q[3];
sx q[3];
rz(1.2866502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6764549) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(-1.6949534) q[2];
rz(2.1220574) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4296221) q[0];
sx q[0];
rz(-1.1638887) q[0];
sx q[0];
rz(0.0042560552) q[0];
rz(2.440522) q[1];
sx q[1];
rz(-0.509976) q[1];
sx q[1];
rz(-3.1140936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6228559) q[0];
sx q[0];
rz(-2.7122926) q[0];
sx q[0];
rz(1.1909498) q[0];
x q[1];
rz(0.84553366) q[2];
sx q[2];
rz(-1.4730244) q[2];
sx q[2];
rz(-0.95319437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78464064) q[1];
sx q[1];
rz(-0.36423794) q[1];
sx q[1];
rz(2.0924773) q[1];
x q[2];
rz(1.2345838) q[3];
sx q[3];
rz(-1.5576406) q[3];
sx q[3];
rz(2.1845078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6444401) q[2];
sx q[2];
rz(-1.7915598) q[2];
sx q[2];
rz(0.39786097) q[2];
rz(-2.1286987) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(-0.40867543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48552805) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(-0.84554607) q[0];
rz(-0.0088508765) q[1];
sx q[1];
rz(-2.2933941) q[1];
sx q[1];
rz(-1.1525851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8542957) q[0];
sx q[0];
rz(-2.3919772) q[0];
sx q[0];
rz(-0.95486705) q[0];
rz(-pi) q[1];
rz(2.8834226) q[2];
sx q[2];
rz(-1.1328146) q[2];
sx q[2];
rz(0.5817619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2611672) q[1];
sx q[1];
rz(-2.5486055) q[1];
sx q[1];
rz(-3.0914607) q[1];
rz(0.20219965) q[3];
sx q[3];
rz(-0.31980896) q[3];
sx q[3];
rz(2.5068738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3919966) q[2];
sx q[2];
rz(-2.1139202) q[2];
sx q[2];
rz(1.2006302) q[2];
rz(-2.6642753) q[3];
sx q[3];
rz(-0.47003191) q[3];
sx q[3];
rz(-2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2068127) q[0];
sx q[0];
rz(-0.88345695) q[0];
sx q[0];
rz(-0.90428895) q[0];
rz(-0.9710871) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(-0.28360525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2663735) q[0];
sx q[0];
rz(-3.0689256) q[0];
sx q[0];
rz(-0.22180264) q[0];
x q[1];
rz(-0.43111151) q[2];
sx q[2];
rz(-1.1656802) q[2];
sx q[2];
rz(-0.27027915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1329258) q[1];
sx q[1];
rz(-0.90331339) q[1];
sx q[1];
rz(2.8776566) q[1];
rz(-pi) q[2];
x q[2];
rz(2.141385) q[3];
sx q[3];
rz(-2.3688487) q[3];
sx q[3];
rz(-2.7955713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7859555) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(-0.23814417) q[2];
rz(-0.98020482) q[3];
sx q[3];
rz(-1.9832059) q[3];
sx q[3];
rz(0.75077209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439483) q[0];
sx q[0];
rz(-1.605796) q[0];
sx q[0];
rz(0.407298) q[0];
rz(0.40335718) q[1];
sx q[1];
rz(-2.401001) q[1];
sx q[1];
rz(-2.2743646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969324) q[0];
sx q[0];
rz(-1.6303902) q[0];
sx q[0];
rz(-1.0048578) q[0];
rz(-pi) q[1];
rz(-1.3620141) q[2];
sx q[2];
rz(-2.1174413) q[2];
sx q[2];
rz(-2.7289313) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2867409) q[1];
sx q[1];
rz(-0.84567243) q[1];
sx q[1];
rz(1.9560567) q[1];
rz(-pi) q[2];
rz(3.0147897) q[3];
sx q[3];
rz(-1.276166) q[3];
sx q[3];
rz(1.0607127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0200218) q[2];
sx q[2];
rz(-0.82841221) q[2];
sx q[2];
rz(1.4309481) q[2];
rz(2.6089148) q[3];
sx q[3];
rz(-2.1055652) q[3];
sx q[3];
rz(-1.7238341) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1581887) q[0];
sx q[0];
rz(-2.991365) q[0];
sx q[0];
rz(1.0839373) q[0];
rz(-3.0896507) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(-0.025024978) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48064104) q[0];
sx q[0];
rz(-1.5990851) q[0];
sx q[0];
rz(-2.6599075) q[0];
x q[1];
rz(2.04966) q[2];
sx q[2];
rz(-2.3298878) q[2];
sx q[2];
rz(1.020249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6083303) q[1];
sx q[1];
rz(-0.46583781) q[1];
sx q[1];
rz(1.9598258) q[1];
rz(1.4281565) q[3];
sx q[3];
rz(-1.9933007) q[3];
sx q[3];
rz(1.6137992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0742801) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(2.32453) q[2];
rz(-0.95064154) q[3];
sx q[3];
rz(-2.1248415) q[3];
sx q[3];
rz(-1.1869717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.10512146) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(0.60086077) q[0];
rz(0.034424456) q[1];
sx q[1];
rz(-1.333678) q[1];
sx q[1];
rz(1.5404125) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8668233) q[0];
sx q[0];
rz(-2.1708198) q[0];
sx q[0];
rz(-2.4680544) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95048423) q[2];
sx q[2];
rz(-2.1455539) q[2];
sx q[2];
rz(1.7983914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.850463) q[1];
sx q[1];
rz(-0.4402059) q[1];
sx q[1];
rz(-3.090211) q[1];
rz(-pi) q[2];
rz(0.65798379) q[3];
sx q[3];
rz(-2.4870928) q[3];
sx q[3];
rz(2.2006048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2711266) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(1.6477443) q[2];
rz(-0.77914023) q[3];
sx q[3];
rz(-1.6611764) q[3];
sx q[3];
rz(-0.65472764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036024) q[0];
sx q[0];
rz(-1.0810738) q[0];
sx q[0];
rz(2.3533452) q[0];
rz(1.8986656) q[1];
sx q[1];
rz(-1.582076) q[1];
sx q[1];
rz(2.8020614) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57836048) q[0];
sx q[0];
rz(-1.6867016) q[0];
sx q[0];
rz(1.3770335) q[0];
rz(1.2844769) q[2];
sx q[2];
rz(-1.0272046) q[2];
sx q[2];
rz(1.7329777) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.06188678) q[1];
sx q[1];
rz(-1.2593653) q[1];
sx q[1];
rz(-2.7175886) q[1];
x q[2];
rz(-2.3655534) q[3];
sx q[3];
rz(-2.0668012) q[3];
sx q[3];
rz(-0.96183323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1303611) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.8846903) q[2];
rz(1.0907178) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(-2.2244661) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85482368) q[0];
sx q[0];
rz(-1.4952861) q[0];
sx q[0];
rz(-2.0538034) q[0];
rz(0.65025672) q[1];
sx q[1];
rz(-1.4573263) q[1];
sx q[1];
rz(3.1204209) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7268459) q[0];
sx q[0];
rz(-1.124589) q[0];
sx q[0];
rz(1.3394974) q[0];
rz(-pi) q[1];
rz(1.7672054) q[2];
sx q[2];
rz(-1.1078664) q[2];
sx q[2];
rz(-0.48505515) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0291527) q[1];
sx q[1];
rz(-2.314056) q[1];
sx q[1];
rz(2.854101) q[1];
x q[2];
rz(1.0383426) q[3];
sx q[3];
rz(-1.5335011) q[3];
sx q[3];
rz(-1.9017912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1374958) q[2];
sx q[2];
rz(-1.3429514) q[2];
sx q[2];
rz(-0.98319483) q[2];
rz(2.148597) q[3];
sx q[3];
rz(-1.1119305) q[3];
sx q[3];
rz(-2.907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642333) q[0];
sx q[0];
rz(-2.3491884) q[0];
sx q[0];
rz(1.0607251) q[0];
rz(1.1119153) q[1];
sx q[1];
rz(-0.30525515) q[1];
sx q[1];
rz(0.12259604) q[1];
rz(-1.1373043) q[2];
sx q[2];
rz(-0.52327427) q[2];
sx q[2];
rz(0.56868194) q[2];
rz(-1.5817011) q[3];
sx q[3];
rz(-1.047871) q[3];
sx q[3];
rz(-2.0497608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
