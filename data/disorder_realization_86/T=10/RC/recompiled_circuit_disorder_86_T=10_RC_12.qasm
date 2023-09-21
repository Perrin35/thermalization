OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(-1.4821948) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9528113) q[0];
sx q[0];
rz(-0.45438284) q[0];
sx q[0];
rz(-2.7812468) q[0];
x q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-0.28684068) q[2];
sx q[2];
rz(-1.4925721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(-0.4835101) q[1];
x q[2];
rz(-0.10098884) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(-1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(2.2757754) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022097691) q[0];
sx q[0];
rz(-2.1848328) q[0];
sx q[0];
rz(0.0028145785) q[0];
rz(-pi) q[1];
rz(-0.33032592) q[2];
sx q[2];
rz(-1.0268372) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7746437) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(0.36555396) q[1];
x q[2];
rz(-2.1727932) q[3];
sx q[3];
rz(-2.7235944) q[3];
sx q[3];
rz(0.19817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(2.1988595) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(0.79743687) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(-0.55999666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277498) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(-1.9613128) q[0];
rz(-pi) q[1];
rz(0.067990818) q[2];
sx q[2];
rz(-2.8376841) q[2];
sx q[2];
rz(1.4246724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3126038) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(2.7752084) q[1];
rz(-pi) q[2];
rz(2.0542164) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(-2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-2.9690572) q[2];
rz(-2.1595188) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(2.8570783) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.8428615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38628681) q[0];
sx q[0];
rz(-2.5884429) q[0];
sx q[0];
rz(-2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80438517) q[2];
sx q[2];
rz(-1.569869) q[2];
sx q[2];
rz(-1.550012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97949308) q[1];
sx q[1];
rz(-2.5433308) q[1];
sx q[1];
rz(-1.23566) q[1];
rz(-pi) q[2];
rz(-1.2508568) q[3];
sx q[3];
rz(-0.3443998) q[3];
sx q[3];
rz(-0.26564769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.3866562) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-0.2968266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083858) q[0];
sx q[0];
rz(-2.1413681) q[0];
sx q[0];
rz(-1.5242566) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74283959) q[2];
sx q[2];
rz(-1.4577216) q[2];
sx q[2];
rz(-2.6807221) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2424803) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(-2.2805023) q[1];
rz(1.4818707) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(-0.54642788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(0.60633916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5099112) q[0];
sx q[0];
rz(-1.5874377) q[0];
sx q[0];
rz(0.04583866) q[0];
x q[1];
rz(2.2491127) q[2];
sx q[2];
rz(-1.9024897) q[2];
sx q[2];
rz(-1.734037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(1.9028266) q[1];
rz(-pi) q[2];
rz(1.243152) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(-0.71654191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(2.6334921) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7471874) q[0];
sx q[0];
rz(-2.1571775) q[0];
sx q[0];
rz(-1.2140843) q[0];
rz(-pi) q[1];
rz(0.21364084) q[2];
sx q[2];
rz(-1.5822615) q[2];
sx q[2];
rz(-1.849888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3527457) q[1];
sx q[1];
rz(-2.6585796) q[1];
sx q[1];
rz(0.61139815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1931476) q[3];
sx q[3];
rz(-0.95818633) q[3];
sx q[3];
rz(2.2539504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(-1.0472263) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(1.3887127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1202639) q[0];
sx q[0];
rz(-1.5409894) q[0];
sx q[0];
rz(-2.1304312) q[0];
rz(1.9954761) q[2];
sx q[2];
rz(-2.0600024) q[2];
sx q[2];
rz(2.3236772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3125004) q[1];
sx q[1];
rz(-1.478985) q[1];
sx q[1];
rz(1.5822259) q[1];
x q[2];
rz(-0.29626131) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(2.5833599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6728014) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(-1.1084491) q[0];
rz(1.6782645) q[2];
sx q[2];
rz(-1.6245914) q[2];
sx q[2];
rz(-1.6966284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3840752) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(-2.3149895) q[1];
rz(2.6551412) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1062746) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(-1.4950698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6388355) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(2.9266561) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(2.4187947) q[1];
rz(-pi) q[2];
rz(0.75533112) q[3];
sx q[3];
rz(-0.33077251) q[3];
sx q[3];
rz(1.8473234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-0.4942016) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(-2.0555029) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.833563) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(-2.3446463) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
