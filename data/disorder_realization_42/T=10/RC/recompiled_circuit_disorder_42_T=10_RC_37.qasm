OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(-2.1142168) q[0];
sx q[0];
rz(2.7789814) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61603123) q[0];
sx q[0];
rz(-1.4481059) q[0];
sx q[0];
rz(0.37567715) q[0];
x q[1];
rz(-1.2149493) q[2];
sx q[2];
rz(-1.0094182) q[2];
sx q[2];
rz(1.9422308) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4075549) q[1];
sx q[1];
rz(-1.2734379) q[1];
sx q[1];
rz(-2.107891) q[1];
rz(-1.9699691) q[3];
sx q[3];
rz(-1.7776383) q[3];
sx q[3];
rz(-0.18404657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.1454426) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(2.4157445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1918068) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(1.1508862) q[0];
rz(-0.59402324) q[2];
sx q[2];
rz(-1.960264) q[2];
sx q[2];
rz(0.71690744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-3.1152578) q[1];
rz(-2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(0.089240616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48250616) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008464) q[0];
sx q[0];
rz(-1.8119537) q[0];
sx q[0];
rz(1.7667889) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4586787) q[2];
sx q[2];
rz(-0.80692601) q[2];
sx q[2];
rz(0.57070953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14721522) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(-1.6836402) q[1];
rz(-pi) q[2];
rz(1.521103) q[3];
sx q[3];
rz(-2.4274785) q[3];
sx q[3];
rz(2.7793838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(2.1543489) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(0.12606829) q[0];
rz(-0.18445045) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(1.1674081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8465189) q[0];
sx q[0];
rz(-2.3210038) q[0];
sx q[0];
rz(-1.2942765) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5037751) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(0.7047082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5645204) q[1];
sx q[1];
rz(-2.2153691) q[1];
sx q[1];
rz(1.9034027) q[1];
rz(-0.79951841) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(-2.094401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(2.5100822) q[2];
rz(-2.946092) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1642078) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(0.50267977) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-2.6745093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750077) q[0];
sx q[0];
rz(-1.9420615) q[0];
sx q[0];
rz(-2.9270372) q[0];
rz(-pi) q[1];
rz(-1.3181114) q[2];
sx q[2];
rz(-2.1420896) q[2];
sx q[2];
rz(-2.2342482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6131763) q[1];
sx q[1];
rz(-1.3532552) q[1];
sx q[1];
rz(2.0030641) q[1];
rz(-2.2561982) q[3];
sx q[3];
rz(-1.4476895) q[3];
sx q[3];
rz(-0.76243329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(1.5658437) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.6451947) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841227) q[0];
sx q[0];
rz(-2.5792482) q[0];
sx q[0];
rz(1.9447295) q[0];
rz(-1.1804579) q[2];
sx q[2];
rz(-2.5679553) q[2];
sx q[2];
rz(2.8611285) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.306957) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(1.6389695) q[1];
rz(-1.5835285) q[3];
sx q[3];
rz(-0.85789645) q[3];
sx q[3];
rz(-1.7604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4795586) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.7842402) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-3.1299652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8893338) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(-0.61347368) q[0];
rz(-0.02545698) q[2];
sx q[2];
rz(-1.5426794) q[2];
sx q[2];
rz(0.97719976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7736518) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(1.6163338) q[1];
rz(-pi) q[2];
rz(2.743268) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(2.6834965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(1.9912432) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(-2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(1.4935965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.134682) q[0];
sx q[0];
rz(-1.784076) q[0];
sx q[0];
rz(-0.30506217) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.059842589) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(3.0795385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4920602) q[1];
sx q[1];
rz(-0.82693716) q[1];
sx q[1];
rz(-2.5745113) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0687374) q[3];
sx q[3];
rz(-1.4668674) q[3];
sx q[3];
rz(-2.7285189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(-1.9624306) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0807063) q[0];
sx q[0];
rz(-1.0651677) q[0];
sx q[0];
rz(1.9067184) q[0];
rz(2.6369008) q[2];
sx q[2];
rz(-1.0070649) q[2];
sx q[2];
rz(-1.6187514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0563911) q[1];
sx q[1];
rz(-1.8279652) q[1];
sx q[1];
rz(0.53253865) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37971244) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(0.1212561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(-1.3423086) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(-0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5060447) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(0.062775469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18386545) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(1.136214) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8809001) q[1];
sx q[1];
rz(-1.1076733) q[1];
sx q[1];
rz(-0.63219627) q[1];
rz(-1.2286931) q[3];
sx q[3];
rz(-1.3918687) q[3];
sx q[3];
rz(-0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(-0.87674117) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(-1.0629874) q[2];
sx q[2];
rz(-0.72401902) q[2];
sx q[2];
rz(0.25821092) q[2];
rz(-3.1003351) q[3];
sx q[3];
rz(-1.2357124) q[3];
sx q[3];
rz(0.99832051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];