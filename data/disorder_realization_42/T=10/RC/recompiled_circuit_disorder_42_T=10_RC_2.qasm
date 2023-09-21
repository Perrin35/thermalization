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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5255614) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(2.7659155) q[0];
rz(-pi) q[1];
rz(0.50589675) q[2];
sx q[2];
rz(-2.4873173) q[2];
sx q[2];
rz(-1.8088532) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33597782) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(-0.34259818) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9176641) q[3];
sx q[3];
rz(-1.9609946) q[3];
sx q[3];
rz(1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.99615) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23913357) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(-0.59160561) q[0];
x q[1];
rz(0.59402324) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(0.71690744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4821856) q[1];
sx q[1];
rz(-0.30989753) q[1];
sx q[1];
rz(-1.6531497) q[1];
x q[2];
rz(0.62249448) q[3];
sx q[3];
rz(-1.3982292) q[3];
sx q[3];
rz(-1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(2.144311) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230351) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(-1.1285271) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(0.9054786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14074621) q[0];
sx q[0];
rz(-1.8119537) q[0];
sx q[0];
rz(1.7667889) q[0];
rz(2.3891874) q[2];
sx q[2];
rz(-1.89626) q[2];
sx q[2];
rz(0.67093713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9943774) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(1.4579525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6204897) q[3];
sx q[3];
rz(-0.71411413) q[3];
sx q[3];
rz(2.7793838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(0.9872438) q[2];
rz(-3.0958214) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(1.1674081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0569939) q[0];
sx q[0];
rz(-1.7718678) q[0];
sx q[0];
rz(0.76954557) q[0];
rz(-pi) q[1];
rz(-2.8235769) q[2];
sx q[2];
rz(-2.615775) q[2];
sx q[2];
rz(1.9981245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(-2.4697484) q[1];
rz(0.64658029) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(-0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(2.946092) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53411667) q[0];
sx q[0];
rz(-0.42629888) q[0];
sx q[0];
rz(-1.0703474) q[0];
rz(-pi) q[1];
rz(-1.3181114) q[2];
sx q[2];
rz(-0.99950302) q[2];
sx q[2];
rz(2.2342482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0568697) q[1];
sx q[1];
rz(-1.9922171) q[1];
sx q[1];
rz(-0.23878581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88539447) q[3];
sx q[3];
rz(-1.4476895) q[3];
sx q[3];
rz(-2.3791594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(-1.5658437) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841227) q[0];
sx q[0];
rz(-0.5623445) q[0];
sx q[0];
rz(-1.9447295) q[0];
rz(-0.24106579) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(-2.966553) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8346356) q[1];
sx q[1];
rz(-2.7863057) q[1];
sx q[1];
rz(-1.5026232) q[1];
rz(-2.4286527) q[3];
sx q[3];
rz(-1.561165) q[3];
sx q[3];
rz(2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.5007639) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-1.0923882) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(0.011627442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2627633) q[0];
sx q[0];
rz(-1.3063523) q[0];
sx q[0];
rz(0.39475616) q[0];
rz(-pi) q[1];
rz(-1.5989223) q[2];
sx q[2];
rz(-1.5962432) q[2];
sx q[2];
rz(-0.59288073) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9776476) q[1];
sx q[1];
rz(-1.5471336) q[1];
sx q[1];
rz(-2.1169099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2467975) q[3];
sx q[3];
rz(-0.9231336) q[3];
sx q[3];
rz(3.0917633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-0.28798506) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
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
rz(-2.8312254) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(-1.2121375) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(1.4935965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15570116) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(-0.62472384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9940894) q[2];
sx q[2];
rz(-2.9967542) q[2];
sx q[2];
rz(-2.7763979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89477506) q[1];
sx q[1];
rz(-2.2404039) q[1];
sx q[1];
rz(2.0991904) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7842403) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(-1.3444927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.2517694) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060886325) q[0];
sx q[0];
rz(-2.076425) q[0];
sx q[0];
rz(1.9067184) q[0];
rz(-pi) q[1];
rz(2.1962677) q[2];
sx q[2];
rz(-1.9918459) q[2];
sx q[2];
rz(-0.23907121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8930248) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(-0.47793169) q[1];
rz(2.7618802) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(0.1212561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(2.9163196) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(0.061696079) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(-1.0888938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5060447) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(-3.0788172) q[0];
x q[1];
rz(-0.91572362) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(-0.54619782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2606925) q[1];
sx q[1];
rz(-1.1076733) q[1];
sx q[1];
rz(-2.5093964) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9519134) q[3];
sx q[3];
rz(-1.2343725) q[3];
sx q[3];
rz(1.3996901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(-0.87674117) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4927647) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(1.0629874) q[2];
sx q[2];
rz(-2.4175736) q[2];
sx q[2];
rz(-2.8833817) q[2];
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