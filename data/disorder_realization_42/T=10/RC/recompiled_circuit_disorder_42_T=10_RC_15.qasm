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
rz(4.1689685) q[0];
sx q[0];
rz(12.203759) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1385961) q[0];
sx q[0];
rz(-1.1980822) q[0];
sx q[0];
rz(1.4390104) q[0];
rz(-pi) q[1];
rz(-1.2149493) q[2];
sx q[2];
rz(-2.1321745) q[2];
sx q[2];
rz(-1.9422308) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33597782) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(-2.7989945) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9699691) q[3];
sx q[3];
rz(-1.3639543) q[3];
sx q[3];
rz(-2.9575461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(0.72584814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23913357) q[0];
sx q[0];
rz(-0.67507889) q[0];
sx q[0];
rz(0.59160561) q[0];
rz(-pi) q[1];
rz(2.0306573) q[2];
sx q[2];
rz(-1.0265372) q[2];
sx q[2];
rz(-2.0366675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.010163807) q[1];
sx q[1];
rz(-1.5958852) q[1];
sx q[1];
rz(-1.8797092) q[1];
x q[2];
rz(-2.5190982) q[3];
sx q[3];
rz(-1.3982292) q[3];
sx q[3];
rz(-1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(0.99728161) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230351) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(2.2361141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008464) q[0];
sx q[0];
rz(-1.329639) q[0];
sx q[0];
rz(-1.3748037) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3891874) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(-2.4706555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29618759) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(0.097464949) q[1];
x q[2];
rz(-1.6204897) q[3];
sx q[3];
rz(-0.71411413) q[3];
sx q[3];
rz(0.36220887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(-1.9741845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4520893) q[0];
sx q[0];
rz(-0.79012442) q[0];
sx q[0];
rz(2.8566314) q[0];
rz(-pi) q[1];
rz(-1.7502968) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(-2.3617982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0988401) q[1];
sx q[1];
rz(-0.71430695) q[1];
sx q[1];
rz(-2.7318098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4950124) q[3];
sx q[3];
rz(-2.1255891) q[3];
sx q[3];
rz(-0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(2.5100822) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-0.41803944) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(0.46708333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53411667) q[0];
sx q[0];
rz(-0.42629888) q[0];
sx q[0];
rz(1.0703474) q[0];
x q[1];
rz(2.7706625) q[2];
sx q[2];
rz(-0.61894722) q[2];
sx q[2];
rz(-1.3528454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.621532) q[1];
sx q[1];
rz(-0.48081765) q[1];
sx q[1];
rz(2.0562999) q[1];
rz(-1.763837) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(0.65929268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6749394) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(-0.28218937) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(0.085993275) q[0];
rz(-1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.6451947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2917568) q[0];
sx q[0];
rz(-2.09016) q[0];
sx q[0];
rz(-0.22626466) q[0];
x q[1];
rz(-1.9611347) q[2];
sx q[2];
rz(-2.5679553) q[2];
sx q[2];
rz(-2.8611285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9073346) q[1];
sx q[1];
rz(-1.9252216) q[1];
sx q[1];
rz(0.025269421) q[1];
x q[2];
rz(2.4286527) q[3];
sx q[3];
rz(-1.5804277) q[3];
sx q[3];
rz(2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(3.1299652) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893338) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(2.528119) q[0];
x q[1];
rz(-0.83519148) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(1.4284301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7736518) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(-1.6163338) q[1];
rz(-pi) q[2];
rz(-1.8947951) q[3];
sx q[3];
rz(-0.9231336) q[3];
sx q[3];
rz(3.0917633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0069107) q[0];
sx q[0];
rz(-1.3575166) q[0];
sx q[0];
rz(-0.30506217) q[0];
rz(-pi) q[1];
rz(3.0817501) q[2];
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
rz(-pi/2) q[0];
rz(-0.89477506) q[1];
sx q[1];
rz(-0.90118876) q[1];
sx q[1];
rz(2.0991904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11843189) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-1.1008319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.9624306) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.8633435) q[0];
sx q[0];
rz(0.53036687) q[0];
x q[1];
rz(2.6369008) q[2];
sx q[2];
rz(-2.1345277) q[2];
sx q[2];
rz(1.6187514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8930248) q[1];
sx q[1];
rz(-0.5859468) q[1];
sx q[1];
rz(-0.47793169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9411373) q[3];
sx q[3];
rz(-2.5459873) q[3];
sx q[3];
rz(-0.60048088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(2.1098095) q[2];
rz(-1.3423086) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(-1.306698) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93931224) q[0];
sx q[0];
rz(-1.6334403) q[0];
sx q[0];
rz(-1.63562) q[0];
rz(-0.91572362) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(2.5953948) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(2.1250493) q[1];
rz(1.0763361) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6488279) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(0.4060612) q[2];
sx q[2];
rz(-2.1880697) q[2];
sx q[2];
rz(-0.38068117) q[2];
rz(-0.041257507) q[3];
sx q[3];
rz(-1.9058803) q[3];
sx q[3];
rz(-2.1432721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];