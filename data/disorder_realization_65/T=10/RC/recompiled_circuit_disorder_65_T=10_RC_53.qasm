OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(-2.8595964) q[0];
x q[1];
rz(-1.7311814) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(-0.10057848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97779146) q[1];
sx q[1];
rz(-1.4492387) q[1];
sx q[1];
rz(-1.4896643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53332897) q[3];
sx q[3];
rz(-1.9682069) q[3];
sx q[3];
rz(-2.1189789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.6404023) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06022913) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(1.5958022) q[0];
rz(-1.00374) q[2];
sx q[2];
rz(-0.54076414) q[2];
sx q[2];
rz(-0.31976779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3850296) q[1];
sx q[1];
rz(-2.2379413) q[1];
sx q[1];
rz(1.7900311) q[1];
rz(-pi) q[2];
rz(0.0060175671) q[3];
sx q[3];
rz(-1.179751) q[3];
sx q[3];
rz(-2.1276377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(2.2272026) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(-0.93018053) q[0];
x q[1];
rz(0.47567993) q[2];
sx q[2];
rz(-1.0389581) q[2];
sx q[2];
rz(1.735641) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55206628) q[1];
sx q[1];
rz(-0.30391903) q[1];
sx q[1];
rz(2.9986831) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91434874) q[3];
sx q[3];
rz(-1.6748866) q[3];
sx q[3];
rz(-2.5641233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(-2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(1.9925041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.944753) q[0];
sx q[0];
rz(-1.6819994) q[0];
sx q[0];
rz(-2.8993594) q[0];
x q[1];
rz(-2.1707702) q[2];
sx q[2];
rz(-0.68471013) q[2];
sx q[2];
rz(0.1240571) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2266352) q[1];
sx q[1];
rz(-1.8640222) q[1];
sx q[1];
rz(-2.3073767) q[1];
rz(-0.37255128) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(-2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(-1.0523798) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28750557) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(2.322312) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7344633) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(1.7083573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21806949) q[1];
sx q[1];
rz(-1.7150208) q[1];
sx q[1];
rz(1.5037392) q[1];
rz(2.103881) q[3];
sx q[3];
rz(-1.8928877) q[3];
sx q[3];
rz(1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(2.172519) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7965308) q[0];
sx q[0];
rz(-0.11632761) q[0];
sx q[0];
rz(-1.2141114) q[0];
rz(-0.078209608) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(-0.90925928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8816996) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(-pi) q[2];
x q[2];
rz(-3*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.4356027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5098269) q[0];
sx q[0];
rz(-1.0477715) q[0];
sx q[0];
rz(2.564389) q[0];
x q[1];
rz(-1.5510467) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(-1.8576647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83511931) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(1.9496099) q[1];
x q[2];
rz(-2.7123659) q[3];
sx q[3];
rz(-1.1042522) q[3];
sx q[3];
rz(-2.6847117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4009565) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.7763604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505969) q[0];
sx q[0];
rz(-0.15582514) q[0];
sx q[0];
rz(-2.0777006) q[0];
rz(-pi) q[1];
rz(-3.119602) q[2];
sx q[2];
rz(-1.4430265) q[2];
sx q[2];
rz(0.86554229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(-0.53175064) q[1];
rz(1.385231) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5101178) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(2.5328013) q[0];
x q[1];
rz(-1.5827492) q[2];
sx q[2];
rz(-2.0593615) q[2];
sx q[2];
rz(1.919463) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9388705) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(2.2321781) q[1];
rz(-pi) q[2];
rz(-0.68848916) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(1.1779932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(2.6031301) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.753153) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.569012) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(1.5241745) q[0];
x q[1];
rz(1.795904) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(3.0181146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.888615) q[1];
sx q[1];
rz(-1.9441609) q[1];
sx q[1];
rz(3.0967767) q[1];
rz(-0.84382236) q[3];
sx q[3];
rz(-2.3355964) q[3];
sx q[3];
rz(2.3674915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(-0.11001982) q[2];
sx q[2];
rz(-1.5129871) q[2];
sx q[2];
rz(-1.5469698) q[2];
rz(0.31717832) q[3];
sx q[3];
rz(-2.2472897) q[3];
sx q[3];
rz(-0.030285611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
