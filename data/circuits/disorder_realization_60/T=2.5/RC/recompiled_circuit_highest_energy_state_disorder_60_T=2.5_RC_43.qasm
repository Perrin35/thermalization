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
rz(0.62975878) q[0];
sx q[0];
rz(3.44343) q[0];
sx q[0];
rz(11.259196) q[0];
rz(-0.52840003) q[1];
sx q[1];
rz(-0.83146787) q[1];
sx q[1];
rz(0.45933476) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6469) q[0];
sx q[0];
rz(-1.6199058) q[0];
sx q[0];
rz(1.8270593) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1856543) q[2];
sx q[2];
rz(-2.8070076) q[2];
sx q[2];
rz(2.0841887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85831538) q[1];
sx q[1];
rz(-2.1766571) q[1];
sx q[1];
rz(1.359904) q[1];
rz(1.5422056) q[3];
sx q[3];
rz(-2.4725879) q[3];
sx q[3];
rz(-2.1137848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.022973148) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(0.55707923) q[2];
rz(2.6856375) q[3];
sx q[3];
rz(-0.7619226) q[3];
sx q[3];
rz(0.60321155) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93422455) q[0];
sx q[0];
rz(-0.71894431) q[0];
sx q[0];
rz(-1.7508605) q[0];
rz(-1.3181744) q[1];
sx q[1];
rz(-1.7134106) q[1];
sx q[1];
rz(-0.21600977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586313) q[0];
sx q[0];
rz(-0.68136946) q[0];
sx q[0];
rz(1.0555223) q[0];
rz(-pi) q[1];
rz(-1.6902709) q[2];
sx q[2];
rz(-1.9873957) q[2];
sx q[2];
rz(2.4753776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4281126) q[1];
sx q[1];
rz(-2.0779582) q[1];
sx q[1];
rz(-0.79913862) q[1];
rz(-1.2162105) q[3];
sx q[3];
rz(-0.56963339) q[3];
sx q[3];
rz(2.5550752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4231437) q[2];
sx q[2];
rz(-2.7687912) q[2];
sx q[2];
rz(-3.0917523) q[2];
rz(-0.54346624) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922888) q[0];
sx q[0];
rz(-1.1019305) q[0];
sx q[0];
rz(-2.1732543) q[0];
rz(-3.1210506) q[1];
sx q[1];
rz(-1.3803218) q[1];
sx q[1];
rz(-1.7113908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50043369) q[0];
sx q[0];
rz(-0.64167385) q[0];
sx q[0];
rz(-0.037317567) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2015351) q[2];
sx q[2];
rz(-2.103269) q[2];
sx q[2];
rz(0.42064253) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9544977) q[1];
sx q[1];
rz(-1.375838) q[1];
sx q[1];
rz(-2.7379237) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3768206) q[3];
sx q[3];
rz(-1.6018036) q[3];
sx q[3];
rz(3.0135807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7057544) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(1.7819116) q[2];
rz(-2.6333574) q[3];
sx q[3];
rz(-1.5225007) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085676) q[0];
sx q[0];
rz(-0.29656947) q[0];
sx q[0];
rz(-1.0067518) q[0];
rz(-0.83167568) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(-1.1078018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9356954) q[0];
sx q[0];
rz(-0.95608053) q[0];
sx q[0];
rz(0.45386916) q[0];
rz(-pi) q[1];
rz(-0.37563373) q[2];
sx q[2];
rz(-1.3357329) q[2];
sx q[2];
rz(-1.7441074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2017758) q[1];
sx q[1];
rz(-0.30758938) q[1];
sx q[1];
rz(-1.5053476) q[1];
rz(-1.8516225) q[3];
sx q[3];
rz(-1.1410332) q[3];
sx q[3];
rz(2.6971779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45336777) q[2];
sx q[2];
rz(-2.0757165) q[2];
sx q[2];
rz(-0.31309703) q[2];
rz(-0.39737663) q[3];
sx q[3];
rz(-1.4704967) q[3];
sx q[3];
rz(0.1853005) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6322286) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(-0.46982345) q[1];
sx q[1];
rz(-1.3146105) q[1];
sx q[1];
rz(-2.125461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7604664) q[0];
sx q[0];
rz(-2.010049) q[0];
sx q[0];
rz(0.45760568) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26922981) q[2];
sx q[2];
rz(-1.2310264) q[2];
sx q[2];
rz(-1.4087848) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.092246902) q[1];
sx q[1];
rz(-1.3692442) q[1];
sx q[1];
rz(-1.818403) q[1];
x q[2];
rz(-1.4757094) q[3];
sx q[3];
rz(-0.32470266) q[3];
sx q[3];
rz(-1.8717757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4887345) q[2];
sx q[2];
rz(-0.32565871) q[2];
sx q[2];
rz(-0.48750901) q[2];
rz(2.2440535) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(2.0770309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(0.94051903) q[0];
sx q[0];
rz(-2.2036393) q[0];
sx q[0];
rz(2.4679389) q[0];
rz(1.9081839) q[1];
sx q[1];
rz(-1.0181095) q[1];
sx q[1];
rz(0.76603755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097283) q[0];
sx q[0];
rz(-1.0610631) q[0];
sx q[0];
rz(-2.8317189) q[0];
x q[1];
rz(0.17579097) q[2];
sx q[2];
rz(-1.5002709) q[2];
sx q[2];
rz(-0.27276892) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8173816) q[1];
sx q[1];
rz(-1.7158475) q[1];
sx q[1];
rz(1.0853898) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30575846) q[3];
sx q[3];
rz(-1.7019203) q[3];
sx q[3];
rz(2.6763889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.172714) q[2];
sx q[2];
rz(-1.6438899) q[2];
sx q[2];
rz(2.3806351) q[2];
rz(-2.739665) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(-0.5886122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589979) q[0];
sx q[0];
rz(-0.75356475) q[0];
sx q[0];
rz(1.6931417) q[0];
rz(-0.50085577) q[1];
sx q[1];
rz(-1.8468937) q[1];
sx q[1];
rz(-2.3282611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091917) q[0];
sx q[0];
rz(-1.4228357) q[0];
sx q[0];
rz(-2.7514821) q[0];
rz(-pi) q[1];
rz(-0.11946875) q[2];
sx q[2];
rz(-1.0257693) q[2];
sx q[2];
rz(1.8785005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5593181) q[1];
sx q[1];
rz(-0.46119237) q[1];
sx q[1];
rz(2.8282437) q[1];
rz(-pi) q[2];
x q[2];
rz(0.058815033) q[3];
sx q[3];
rz(-1.6928738) q[3];
sx q[3];
rz(2.5466145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8616051) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-0.32312265) q[2];
rz(1.7396287) q[3];
sx q[3];
rz(-1.8596545) q[3];
sx q[3];
rz(-1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0793656) q[0];
sx q[0];
rz(-2.0923738) q[0];
sx q[0];
rz(0.76612377) q[0];
rz(0.8194204) q[1];
sx q[1];
rz(-2.1111646) q[1];
sx q[1];
rz(-1.6105509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9510244) q[0];
sx q[0];
rz(-1.6389209) q[0];
sx q[0];
rz(-0.22916746) q[0];
rz(2.4862637) q[2];
sx q[2];
rz(-1.3957784) q[2];
sx q[2];
rz(-1.1879712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30554015) q[1];
sx q[1];
rz(-2.17872) q[1];
sx q[1];
rz(0.72956438) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29115486) q[3];
sx q[3];
rz(-1.2586354) q[3];
sx q[3];
rz(-0.96033421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7649585) q[2];
sx q[2];
rz(-0.19516334) q[2];
sx q[2];
rz(2.7435319) q[2];
rz(-2.5352488) q[3];
sx q[3];
rz(-1.5667934) q[3];
sx q[3];
rz(-2.2280367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23676087) q[0];
sx q[0];
rz(-0.30192152) q[0];
sx q[0];
rz(-0.52445573) q[0];
rz(-0.66871387) q[1];
sx q[1];
rz(-1.5293744) q[1];
sx q[1];
rz(1.761577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3052914) q[0];
sx q[0];
rz(-2.143465) q[0];
sx q[0];
rz(3.089726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8369806) q[2];
sx q[2];
rz(-2.0292943) q[2];
sx q[2];
rz(3.0625513) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.804033) q[1];
sx q[1];
rz(-1.1566678) q[1];
sx q[1];
rz(-2.0284589) q[1];
rz(0.53219019) q[3];
sx q[3];
rz(-1.2423774) q[3];
sx q[3];
rz(-1.4946973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9809208) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(0.46372947) q[2];
rz(1.7175698) q[3];
sx q[3];
rz(-2.0537036) q[3];
sx q[3];
rz(-0.033528479) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(-2.723519) q[0];
rz(0.75830165) q[1];
sx q[1];
rz(-0.98379358) q[1];
sx q[1];
rz(1.8565348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.304628) q[0];
sx q[0];
rz(-1.6462176) q[0];
sx q[0];
rz(-2.4914431) q[0];
rz(-0.94916085) q[2];
sx q[2];
rz(-1.2864603) q[2];
sx q[2];
rz(-1.9380707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1753208) q[1];
sx q[1];
rz(-1.0283677) q[1];
sx q[1];
rz(1.7955417) q[1];
x q[2];
rz(2.140652) q[3];
sx q[3];
rz(-2.7469198) q[3];
sx q[3];
rz(-0.81854023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16984223) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(2.8954835) q[2];
rz(-1.2717815) q[3];
sx q[3];
rz(-1.566794) q[3];
sx q[3];
rz(-1.7790214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590189) q[0];
sx q[0];
rz(-1.4763426) q[0];
sx q[0];
rz(2.939298) q[0];
rz(-0.62494878) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(-0.18904674) q[2];
sx q[2];
rz(-1.3658267) q[2];
sx q[2];
rz(1.7455802) q[2];
rz(-0.20062867) q[3];
sx q[3];
rz(-1.8798141) q[3];
sx q[3];
rz(0.78165913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
