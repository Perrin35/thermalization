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
rz(-0.15859088) q[0];
sx q[0];
rz(2.6245485) q[0];
sx q[0];
rz(11.256097) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(0.96460834) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91473305) q[0];
sx q[0];
rz(-2.4977685) q[0];
sx q[0];
rz(-2.9152246) q[0];
rz(-pi) q[1];
rz(-2.7959441) q[2];
sx q[2];
rz(-1.7549577) q[2];
sx q[2];
rz(3.057827) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55935153) q[1];
sx q[1];
rz(-2.3326477) q[1];
sx q[1];
rz(-1.9870583) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89972605) q[3];
sx q[3];
rz(-1.4802855) q[3];
sx q[3];
rz(-2.8495726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0059119314) q[2];
sx q[2];
rz(-0.19401208) q[2];
sx q[2];
rz(1.4646336) q[2];
rz(-1.8320734) q[3];
sx q[3];
rz(-2.194761) q[3];
sx q[3];
rz(2.8215698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784073) q[0];
sx q[0];
rz(-2.0018556) q[0];
sx q[0];
rz(1.3673258) q[0];
rz(-1.4890081) q[1];
sx q[1];
rz(-0.79843489) q[1];
sx q[1];
rz(2.8489825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45918834) q[0];
sx q[0];
rz(-1.4727797) q[0];
sx q[0];
rz(0.3353663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9797249) q[2];
sx q[2];
rz(-1.6965995) q[2];
sx q[2];
rz(0.51728546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.039174883) q[1];
sx q[1];
rz(-2.8264649) q[1];
sx q[1];
rz(2.4246895) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5047362) q[3];
sx q[3];
rz(-1.3050029) q[3];
sx q[3];
rz(-0.059826033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9802398) q[2];
sx q[2];
rz(-2.5148401) q[2];
sx q[2];
rz(-0.18388595) q[2];
rz(-2.6889177) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(0.12990738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9420796) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-0.80514258) q[0];
rz(1.0792271) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(1.315518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53213648) q[0];
sx q[0];
rz(-2.2972496) q[0];
sx q[0];
rz(1.491602) q[0];
rz(-pi) q[1];
rz(-1.075885) q[2];
sx q[2];
rz(-2.333771) q[2];
sx q[2];
rz(-3.0095095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6650369) q[1];
sx q[1];
rz(-1.2893234) q[1];
sx q[1];
rz(-0.095515619) q[1];
rz(1.4844378) q[3];
sx q[3];
rz(-1.8724073) q[3];
sx q[3];
rz(-2.3197966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9341854) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(0.54971131) q[2];
rz(3.1304729) q[3];
sx q[3];
rz(-0.79922262) q[3];
sx q[3];
rz(1.3219249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84871197) q[0];
sx q[0];
rz(-2.1826545) q[0];
sx q[0];
rz(2.8380561) q[0];
rz(-2.983298) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(2.1319481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5925668) q[0];
sx q[0];
rz(-2.0011304) q[0];
sx q[0];
rz(0.31220147) q[0];
rz(-pi) q[1];
rz(0.60229566) q[2];
sx q[2];
rz(-1.8297075) q[2];
sx q[2];
rz(0.4073172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6827439) q[1];
sx q[1];
rz(-1.413133) q[1];
sx q[1];
rz(-2.661652) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49780881) q[3];
sx q[3];
rz(-1.5524929) q[3];
sx q[3];
rz(-2.8895484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2046854) q[2];
sx q[2];
rz(-0.82131177) q[2];
sx q[2];
rz(2.8724907) q[2];
rz(1.4540295) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(-1.5627741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8708385) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(-2.7794072) q[0];
rz(-2.4902792) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(-1.6168894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37446132) q[0];
sx q[0];
rz(-0.62333019) q[0];
sx q[0];
rz(0.222739) q[0];
rz(-pi) q[1];
rz(-0.94793041) q[2];
sx q[2];
rz(-1.6375223) q[2];
sx q[2];
rz(0.78363505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73495364) q[1];
sx q[1];
rz(-1.2124774) q[1];
sx q[1];
rz(2.1291158) q[1];
rz(-pi) q[2];
rz(1.175143) q[3];
sx q[3];
rz(-0.86651245) q[3];
sx q[3];
rz(0.3053529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1093971) q[2];
sx q[2];
rz(-2.1869662) q[2];
sx q[2];
rz(1.3738013) q[2];
rz(-2.7023756) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(0.60905987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1010901) q[0];
sx q[0];
rz(-2.4093565) q[0];
sx q[0];
rz(-2.8771583) q[0];
rz(2.4624372) q[1];
sx q[1];
rz(-2.2676088) q[1];
sx q[1];
rz(0.98495475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20709461) q[0];
sx q[0];
rz(-2.6885928) q[0];
sx q[0];
rz(-0.82243246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0206775) q[2];
sx q[2];
rz(-1.3123672) q[2];
sx q[2];
rz(2.8095736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9912274) q[1];
sx q[1];
rz(-0.62365198) q[1];
sx q[1];
rz(-0.6536478) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1285176) q[3];
sx q[3];
rz(-1.7566214) q[3];
sx q[3];
rz(3.1215661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(-0.43231535) q[2];
rz(1.013422) q[3];
sx q[3];
rz(-1.5347967) q[3];
sx q[3];
rz(-2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6559615) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(1.1095169) q[0];
rz(-0.79311496) q[1];
sx q[1];
rz(-1.7879281) q[1];
sx q[1];
rz(-1.5951593) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4594134) q[0];
sx q[0];
rz(-2.2930682) q[0];
sx q[0];
rz(0.79011495) q[0];
x q[1];
rz(2.903669) q[2];
sx q[2];
rz(-0.73582651) q[2];
sx q[2];
rz(3.0032681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82779694) q[1];
sx q[1];
rz(-0.89875752) q[1];
sx q[1];
rz(-2.8845877) q[1];
rz(-pi) q[2];
rz(2.4251806) q[3];
sx q[3];
rz(-1.6509112) q[3];
sx q[3];
rz(0.83601213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82211295) q[2];
sx q[2];
rz(-2.959368) q[2];
sx q[2];
rz(1.7721843) q[2];
rz(-0.93044126) q[3];
sx q[3];
rz(-1.3481827) q[3];
sx q[3];
rz(-1.6377009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9419788) q[0];
sx q[0];
rz(-1.9604585) q[0];
sx q[0];
rz(2.8379295) q[0];
rz(-2.1619201) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(2.7365275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41202085) q[0];
sx q[0];
rz(-1.7639909) q[0];
sx q[0];
rz(-0.3839107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.943739) q[2];
sx q[2];
rz(-2.6952546) q[2];
sx q[2];
rz(2.9624903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0822337) q[1];
sx q[1];
rz(-0.59525604) q[1];
sx q[1];
rz(-0.31945503) q[1];
x q[2];
rz(-1.1478506) q[3];
sx q[3];
rz(-0.4627403) q[3];
sx q[3];
rz(3.1266318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27641174) q[2];
sx q[2];
rz(-2.1750735) q[2];
sx q[2];
rz(0.69918862) q[2];
rz(2.3326323) q[3];
sx q[3];
rz(-1.304108) q[3];
sx q[3];
rz(-2.8185524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.9843822) q[0];
sx q[0];
rz(-0.13741563) q[0];
rz(-3.0925062) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415599) q[0];
sx q[0];
rz(-1.5440436) q[0];
sx q[0];
rz(0.49515611) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069641308) q[2];
sx q[2];
rz(-0.17270522) q[2];
sx q[2];
rz(2.791419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6235828) q[1];
sx q[1];
rz(-1.8933663) q[1];
sx q[1];
rz(-2.0913893) q[1];
rz(-pi) q[2];
rz(-0.089506702) q[3];
sx q[3];
rz(-1.6920768) q[3];
sx q[3];
rz(2.3119777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55234838) q[2];
sx q[2];
rz(-1.2184703) q[2];
sx q[2];
rz(2.153896) q[2];
rz(-0.47932953) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.669303) q[0];
sx q[0];
rz(-2.9030114) q[0];
sx q[0];
rz(-3.0008089) q[0];
rz(-2.7760778) q[1];
sx q[1];
rz(-0.82217685) q[1];
sx q[1];
rz(-0.64868322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012039579) q[0];
sx q[0];
rz(-0.7551935) q[0];
sx q[0];
rz(2.6213403) q[0];
rz(-pi) q[1];
rz(-3.0200483) q[2];
sx q[2];
rz(-0.76833188) q[2];
sx q[2];
rz(0.8476846) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71350559) q[1];
sx q[1];
rz(-2.1660342) q[1];
sx q[1];
rz(1.5028605) q[1];
x q[2];
rz(-1.8126103) q[3];
sx q[3];
rz(-1.2529272) q[3];
sx q[3];
rz(-2.1235025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4623798) q[2];
sx q[2];
rz(-1.0054192) q[2];
sx q[2];
rz(3.0688378) q[2];
rz(0.66271979) q[3];
sx q[3];
rz(-1.8279671) q[3];
sx q[3];
rz(2.976118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.0543906) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(1.5351334) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(-1.6894111) q[2];
sx q[2];
rz(-1.7977503) q[2];
sx q[2];
rz(1.1888421) q[2];
rz(2.2519464) q[3];
sx q[3];
rz(-1.8915265) q[3];
sx q[3];
rz(-1.693207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
