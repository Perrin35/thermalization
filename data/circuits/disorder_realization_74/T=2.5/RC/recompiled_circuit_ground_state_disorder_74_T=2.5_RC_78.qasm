OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7862406) q[0];
sx q[0];
rz(-2.8347637) q[0];
sx q[0];
rz(2.8428349) q[0];
rz(2.4755251) q[1];
sx q[1];
rz(-2.5112285) q[1];
sx q[1];
rz(-1.4730374) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10160343) q[0];
sx q[0];
rz(-0.77058661) q[0];
sx q[0];
rz(-1.5250852) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1236727) q[2];
sx q[2];
rz(-1.8715053) q[2];
sx q[2];
rz(-1.199388) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3949497) q[1];
sx q[1];
rz(-2.267845) q[1];
sx q[1];
rz(0.66956981) q[1];
rz(-pi) q[2];
rz(-1.500559) q[3];
sx q[3];
rz(-1.8459709) q[3];
sx q[3];
rz(0.17632139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.95639688) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(-1.8308651) q[2];
rz(3.0500566) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(-2.6105647) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071534261) q[0];
sx q[0];
rz(-2.0508843) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(1.1473038) q[1];
sx q[1];
rz(-0.6310178) q[1];
sx q[1];
rz(2.4400585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.498256) q[0];
sx q[0];
rz(-2.069888) q[0];
sx q[0];
rz(-0.82478158) q[0];
x q[1];
rz(-2.0364136) q[2];
sx q[2];
rz(-1.1934416) q[2];
sx q[2];
rz(-1.5355664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5614723) q[1];
sx q[1];
rz(-2.0144723) q[1];
sx q[1];
rz(0.89718141) q[1];
x q[2];
rz(-1.8789038) q[3];
sx q[3];
rz(-2.1854221) q[3];
sx q[3];
rz(-2.7648787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7445765) q[2];
sx q[2];
rz(-1.1822367) q[2];
sx q[2];
rz(1.5632632) q[2];
rz(0.27966106) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(-0.40951148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(-1.718234) q[0];
rz(-3.0176945) q[1];
sx q[1];
rz(-0.84404498) q[1];
sx q[1];
rz(1.557225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6829872) q[0];
sx q[0];
rz(-1.8667954) q[0];
sx q[0];
rz(-0.83403765) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36735518) q[2];
sx q[2];
rz(-0.88161385) q[2];
sx q[2];
rz(-0.85509461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2899921) q[1];
sx q[1];
rz(-0.25892913) q[1];
sx q[1];
rz(1.2570791) q[1];
rz(-pi) q[2];
rz(0.73488124) q[3];
sx q[3];
rz(-2.2874444) q[3];
sx q[3];
rz(1.411996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0763756) q[2];
sx q[2];
rz(-0.6535483) q[2];
sx q[2];
rz(0.93488133) q[2];
rz(0.30385083) q[3];
sx q[3];
rz(-2.5287718) q[3];
sx q[3];
rz(-2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(2.8745162) q[0];
rz(-3.0067387) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(2.3042302) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1564604) q[0];
sx q[0];
rz(-1.3044668) q[0];
sx q[0];
rz(1.0087031) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58394152) q[2];
sx q[2];
rz(-2.2211233) q[2];
sx q[2];
rz(-3.1027628) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4761047) q[1];
sx q[1];
rz(-1.3766411) q[1];
sx q[1];
rz(-2.9041415) q[1];
rz(-1.5092974) q[3];
sx q[3];
rz(-1.5240747) q[3];
sx q[3];
rz(-0.72544569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59552938) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(-2.9626633) q[2];
rz(-1.9723802) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(2.9235212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021407481) q[0];
sx q[0];
rz(-0.51161259) q[0];
sx q[0];
rz(-2.2688493) q[0];
rz(1.9841638) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(-3.1246368) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.884753) q[0];
sx q[0];
rz(-1.5829979) q[0];
sx q[0];
rz(-2.4069465) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5418607) q[2];
sx q[2];
rz(-1.5790734) q[2];
sx q[2];
rz(3.099583) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53345799) q[1];
sx q[1];
rz(-2.3124586) q[1];
sx q[1];
rz(1.9749179) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26930974) q[3];
sx q[3];
rz(-1.8519173) q[3];
sx q[3];
rz(-0.65015974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0852647) q[2];
sx q[2];
rz(-1.4339829) q[2];
sx q[2];
rz(2.8223574) q[2];
rz(-2.4625835) q[3];
sx q[3];
rz(-1.346799) q[3];
sx q[3];
rz(2.3398248) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(2.4989682) q[0];
rz(-1.7721666) q[1];
sx q[1];
rz(-2.5182928) q[1];
sx q[1];
rz(0.97698897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8606023) q[0];
sx q[0];
rz(-2.9837768) q[0];
sx q[0];
rz(2.4688323) q[0];
rz(-pi) q[1];
rz(1.0426635) q[2];
sx q[2];
rz(-1.1234049) q[2];
sx q[2];
rz(1.7698947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.381666) q[1];
sx q[1];
rz(-1.5569512) q[1];
sx q[1];
rz(-2.2626876) q[1];
x q[2];
rz(1.6518408) q[3];
sx q[3];
rz(-0.5872927) q[3];
sx q[3];
rz(-1.0451661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6151108) q[2];
sx q[2];
rz(-1.4075764) q[2];
sx q[2];
rz(1.2283481) q[2];
rz(0.86722106) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(-0.76798463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7277471) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(0.17284285) q[0];
rz(0.8051644) q[1];
sx q[1];
rz(-0.54904896) q[1];
sx q[1];
rz(-0.13596143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522536) q[0];
sx q[0];
rz(-0.87119188) q[0];
sx q[0];
rz(-0.1564349) q[0];
rz(-pi) q[1];
rz(3.0375541) q[2];
sx q[2];
rz(-2.316873) q[2];
sx q[2];
rz(-2.4989043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6579305) q[1];
sx q[1];
rz(-1.4560946) q[1];
sx q[1];
rz(1.5516267) q[1];
x q[2];
rz(-1.436266) q[3];
sx q[3];
rz(-1.2438275) q[3];
sx q[3];
rz(3.0643058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92676306) q[2];
sx q[2];
rz(-1.5952933) q[2];
sx q[2];
rz(-2.9570441) q[2];
rz(2.9522225) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9823343) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(0.92639297) q[0];
rz(3.1023846) q[1];
sx q[1];
rz(-2.4640633) q[1];
sx q[1];
rz(2.1047986) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568615) q[0];
sx q[0];
rz(-2.7960092) q[0];
sx q[0];
rz(1.0285818) q[0];
rz(-pi) q[1];
rz(-0.37929566) q[2];
sx q[2];
rz(-0.46078983) q[2];
sx q[2];
rz(2.5658105) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1649303) q[1];
sx q[1];
rz(-1.9769985) q[1];
sx q[1];
rz(-1.3066533) q[1];
rz(-pi) q[2];
rz(-2.7091712) q[3];
sx q[3];
rz(-2.7161971) q[3];
sx q[3];
rz(-2.2063125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89256531) q[2];
sx q[2];
rz(-1.1250863) q[2];
sx q[2];
rz(-0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(2.5911205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(1.5193526) q[0];
rz(1.5937357) q[1];
sx q[1];
rz(-0.82031119) q[1];
sx q[1];
rz(-0.69040745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102745) q[0];
sx q[0];
rz(-1.9209492) q[0];
sx q[0];
rz(1.3331593) q[0];
rz(-pi) q[1];
rz(1.4682653) q[2];
sx q[2];
rz(-2.3674115) q[2];
sx q[2];
rz(-0.098086327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6505594) q[1];
sx q[1];
rz(-2.0692252) q[1];
sx q[1];
rz(0.9935324) q[1];
rz(-pi) q[2];
rz(-0.14491187) q[3];
sx q[3];
rz(-1.0365126) q[3];
sx q[3];
rz(2.2458959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4385628) q[2];
sx q[2];
rz(-0.60882336) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(-2.8841833) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(-1.7820057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064706) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(-3.0560793) q[0];
rz(-0.27587786) q[1];
sx q[1];
rz(-2.52067) q[1];
sx q[1];
rz(-1.1891018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10201926) q[0];
sx q[0];
rz(-1.1652526) q[0];
sx q[0];
rz(-1.1723551) q[0];
rz(-1.5618454) q[2];
sx q[2];
rz(-1.2607961) q[2];
sx q[2];
rz(1.6197038) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.117529) q[1];
sx q[1];
rz(-1.4520565) q[1];
sx q[1];
rz(1.1681144) q[1];
x q[2];
rz(-1.4776049) q[3];
sx q[3];
rz(-1.9177327) q[3];
sx q[3];
rz(3.0938397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(-2.6574668) q[2];
rz(3.0020946) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(-2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.6839266) q[1];
sx q[1];
rz(-1.4977581) q[1];
sx q[1];
rz(-1.3118634) q[1];
rz(1.3670078) q[2];
sx q[2];
rz(-2.2098354) q[2];
sx q[2];
rz(-0.60571203) q[2];
rz(0.75705725) q[3];
sx q[3];
rz(-1.5499877) q[3];
sx q[3];
rz(-3.1201759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
