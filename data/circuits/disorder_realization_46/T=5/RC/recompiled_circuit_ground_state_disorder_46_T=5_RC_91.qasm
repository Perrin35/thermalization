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
rz(5.0862105) q[1];
sx q[1];
rz(12.401019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7900766) q[0];
sx q[0];
rz(-0.72203239) q[0];
sx q[0];
rz(2.0603098) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1520391) q[2];
sx q[2];
rz(-1.1606231) q[2];
sx q[2];
rz(2.9048186) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27439427) q[1];
sx q[1];
rz(-1.9748678) q[1];
sx q[1];
rz(-2.1073125) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88460102) q[3];
sx q[3];
rz(-2.7470318) q[3];
sx q[3];
rz(-1.3742428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87236658) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(-1.1671789) q[2];
rz(0.40927467) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(-2.1854713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-2.9874711) q[0];
sx q[0];
rz(-1.3519721) q[0];
rz(-1.6701472) q[1];
sx q[1];
rz(-0.92266005) q[1];
sx q[1];
rz(-2.2460489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038747) q[0];
sx q[0];
rz(-2.8744937) q[0];
sx q[0];
rz(-0.0077203053) q[0];
rz(-2.6673253) q[2];
sx q[2];
rz(-1.6401498) q[2];
sx q[2];
rz(2.4276395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.99341678) q[1];
sx q[1];
rz(-1.1779146) q[1];
sx q[1];
rz(-0.82156374) q[1];
rz(-pi) q[2];
rz(-0.23625613) q[3];
sx q[3];
rz(-1.9801664) q[3];
sx q[3];
rz(0.3796814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6764549) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(-1.4466393) q[2];
rz(-2.1220574) q[3];
sx q[3];
rz(-1.8229702) q[3];
sx q[3];
rz(-0.29729602) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7119706) q[0];
sx q[0];
rz(-1.977704) q[0];
sx q[0];
rz(3.1373366) q[0];
rz(-2.440522) q[1];
sx q[1];
rz(-0.509976) q[1];
sx q[1];
rz(3.1140936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.741318) q[0];
sx q[0];
rz(-1.4158465) q[0];
sx q[0];
rz(1.9727895) q[0];
rz(-pi) q[1];
x q[1];
rz(2.296059) q[2];
sx q[2];
rz(-1.6685682) q[2];
sx q[2];
rz(-0.95319437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78464064) q[1];
sx q[1];
rz(-0.36423794) q[1];
sx q[1];
rz(1.0491153) q[1];
x q[2];
rz(-0.013935907) q[3];
sx q[3];
rz(-1.2346141) q[3];
sx q[3];
rz(-0.61830904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6444401) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(-2.7437317) q[2];
rz(-1.0128939) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(-2.7329172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6560646) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(-2.2960466) q[0];
rz(3.1327418) q[1];
sx q[1];
rz(-2.2933941) q[1];
sx q[1];
rz(-1.1525851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6229079) q[0];
sx q[0];
rz(-2.1605413) q[0];
sx q[0];
rz(-2.6481762) q[0];
rz(-pi) q[1];
rz(-2.0699224) q[2];
sx q[2];
rz(-2.6374657) q[2];
sx q[2];
rz(-2.0029411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8735476) q[1];
sx q[1];
rz(-1.598804) q[1];
sx q[1];
rz(-2.5491879) q[1];
x q[2];
rz(2.939393) q[3];
sx q[3];
rz(-0.31980896) q[3];
sx q[3];
rz(-2.5068738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3919966) q[2];
sx q[2];
rz(-1.0276724) q[2];
sx q[2];
rz(-1.2006302) q[2];
rz(-0.47731733) q[3];
sx q[3];
rz(-0.47003191) q[3];
sx q[3];
rz(-0.70394713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93477997) q[0];
sx q[0];
rz(-0.88345695) q[0];
sx q[0];
rz(0.90428895) q[0];
rz(2.1705056) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(2.8579874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47434092) q[0];
sx q[0];
rz(-1.5867689) q[0];
sx q[0];
rz(-3.0706997) q[0];
rz(2.3432557) q[2];
sx q[2];
rz(-2.5588648) q[2];
sx q[2];
rz(-0.59205627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7212977) q[1];
sx q[1];
rz(-2.4313213) q[1];
sx q[1];
rz(-1.8904449) q[1];
x q[2];
rz(-2.6568703) q[3];
sx q[3];
rz(-0.94281632) q[3];
sx q[3];
rz(-2.064765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35563719) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(0.23814417) q[2];
rz(2.1613878) q[3];
sx q[3];
rz(-1.1583867) q[3];
sx q[3];
rz(2.3908206) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439483) q[0];
sx q[0];
rz(-1.5357966) q[0];
sx q[0];
rz(-0.407298) q[0];
rz(-2.7382355) q[1];
sx q[1];
rz(-0.74059161) q[1];
sx q[1];
rz(2.2743646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3639528) q[0];
sx q[0];
rz(-1.0059851) q[0];
sx q[0];
rz(0.070567957) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7795785) q[2];
sx q[2];
rz(-1.0241514) q[2];
sx q[2];
rz(-0.41266135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7379511) q[1];
sx q[1];
rz(-0.80437901) q[1];
sx q[1];
rz(-2.7405095) q[1];
rz(-3.0147897) q[3];
sx q[3];
rz(-1.8654267) q[3];
sx q[3];
rz(-2.08088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0200218) q[2];
sx q[2];
rz(-2.3131804) q[2];
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
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.983404) q[0];
sx q[0];
rz(-2.991365) q[0];
sx q[0];
rz(1.0839373) q[0];
rz(-0.051941959) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(0.025024978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6609516) q[0];
sx q[0];
rz(-1.5990851) q[0];
sx q[0];
rz(2.6599075) q[0];
rz(-pi) q[1];
rz(2.3228753) q[2];
sx q[2];
rz(-1.229964) q[2];
sx q[2];
rz(-2.9342295) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6083303) q[1];
sx q[1];
rz(-2.6757548) q[1];
sx q[1];
rz(-1.1817668) q[1];
rz(-pi) q[2];
rz(-0.30625108) q[3];
sx q[3];
rz(-0.44455498) q[3];
sx q[3];
rz(1.2769092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0673125) q[2];
sx q[2];
rz(-1.4558027) q[2];
sx q[2];
rz(2.32453) q[2];
rz(-0.95064154) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(-1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0364712) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(0.60086077) q[0];
rz(3.1071682) q[1];
sx q[1];
rz(-1.333678) q[1];
sx q[1];
rz(-1.5404125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87266028) q[0];
sx q[0];
rz(-1.0300228) q[0];
sx q[0];
rz(0.85178225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.41018) q[2];
sx q[2];
rz(-0.81899511) q[2];
sx q[2];
rz(0.42289823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7936777) q[1];
sx q[1];
rz(-1.131212) q[1];
sx q[1];
rz(-1.5949834) q[1];
rz(-pi) q[2];
rz(-0.54564666) q[3];
sx q[3];
rz(-1.1893404) q[3];
sx q[3];
rz(-1.9617401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8704661) q[2];
sx q[2];
rz(-1.2721964) q[2];
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
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43799022) q[0];
sx q[0];
rz(-1.0810738) q[0];
sx q[0];
rz(-0.78824747) q[0];
rz(1.2429271) q[1];
sx q[1];
rz(-1.5595167) q[1];
sx q[1];
rz(2.8020614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718455) q[0];
sx q[0];
rz(-1.3783499) q[0];
sx q[0];
rz(-0.11809531) q[0];
rz(-pi) q[1];
rz(2.7044136) q[2];
sx q[2];
rz(-2.533982) q[2];
sx q[2];
rz(-0.89113441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7701246) q[1];
sx q[1];
rz(-1.9731908) q[1];
sx q[1];
rz(-1.2312908) q[1];
rz(-pi) q[2];
rz(0.65776627) q[3];
sx q[3];
rz(-0.89221803) q[3];
sx q[3];
rz(-0.15746169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.01123151) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(-1.2569024) q[2];
rz(-1.0907178) q[3];
sx q[3];
rz(-1.9382449) q[3];
sx q[3];
rz(0.91712657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85482368) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(1.0877892) q[0];
rz(-2.4913359) q[1];
sx q[1];
rz(-1.4573263) q[1];
sx q[1];
rz(3.1204209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0567704) q[0];
sx q[0];
rz(-2.6426043) q[0];
sx q[0];
rz(-0.4468687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4707321) q[2];
sx q[2];
rz(-1.3952878) q[2];
sx q[2];
rz(-0.99711768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0291527) q[1];
sx q[1];
rz(-2.314056) q[1];
sx q[1];
rz(-0.2874916) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4974277) q[3];
sx q[3];
rz(-2.6079599) q[3];
sx q[3];
rz(-0.39419202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1374958) q[2];
sx q[2];
rz(-1.3429514) q[2];
sx q[2];
rz(0.98319483) q[2];
rz(-0.99299562) q[3];
sx q[3];
rz(-1.1119305) q[3];
sx q[3];
rz(-2.907311) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642333) q[0];
sx q[0];
rz(-0.79240427) q[0];
sx q[0];
rz(-2.0808676) q[0];
rz(-1.1119153) q[1];
sx q[1];
rz(-2.8363375) q[1];
sx q[1];
rz(-3.0189966) q[1];
rz(1.1373043) q[2];
sx q[2];
rz(-2.6183184) q[2];
sx q[2];
rz(-2.5729107) q[2];
rz(-2.6186416) q[3];
sx q[3];
rz(-1.5613489) q[3];
sx q[3];
rz(-0.47351826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
