OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(-2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0029966) q[0];
sx q[0];
rz(-1.9435104) q[0];
sx q[0];
rz(-1.7025823) q[0];
x q[1];
rz(1.2149493) q[2];
sx q[2];
rz(-1.0094182) q[2];
sx q[2];
rz(1.1993619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8056148) q[1];
sx q[1];
rz(-2.0819547) q[1];
sx q[1];
rz(-2.7989945) q[1];
rz(-pi) q[2];
rz(2.0658675) q[3];
sx q[3];
rz(-2.6945811) q[3];
sx q[3];
rz(-1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.5216924) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(0.72584814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94978588) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(-1.9907065) q[0];
rz(-pi) q[1];
rz(0.63273301) q[2];
sx q[2];
rz(-2.4443812) q[2];
sx q[2];
rz(0.34174191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5729546) q[1];
sx q[1];
rz(-1.2619839) q[1];
sx q[1];
rz(3.1152578) q[1];
rz(-pi) q[2];
rz(-0.62249448) q[3];
sx q[3];
rz(-1.3982292) q[3];
sx q[3];
rz(1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41855758) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(1.1285271) q[0];
rz(1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.758924) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(-0.24567901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3891874) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(2.4706555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29618759) q[1];
sx q[1];
rz(-0.86128174) q[1];
sx q[1];
rz(-0.097464949) q[1];
x q[2];
rz(-3.0985673) q[3];
sx q[3];
rz(-0.85775162) q[3];
sx q[3];
rz(-2.713664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(0.18019095) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4520893) q[0];
sx q[0];
rz(-2.3514682) q[0];
sx q[0];
rz(-2.8566314) q[0];
x q[1];
rz(2.8235769) q[2];
sx q[2];
rz(-0.52581767) q[2];
sx q[2];
rz(1.9981245) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5645204) q[1];
sx q[1];
rz(-0.92622354) q[1];
sx q[1];
rz(1.9034027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91058369) q[3];
sx q[3];
rz(-2.108421) q[3];
sx q[3];
rz(2.0349353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(0.19550066) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(2.6389129) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-2.6745093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53411667) q[0];
sx q[0];
rz(-2.7152938) q[0];
sx q[0];
rz(2.0712453) q[0];
rz(2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(2.6168602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.621532) q[1];
sx q[1];
rz(-2.660775) q[1];
sx q[1];
rz(-2.0562999) q[1];
rz(2.2561982) q[3];
sx q[3];
rz(-1.4476895) q[3];
sx q[3];
rz(-2.3791594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6749394) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-0.085993275) q[0];
rz(-1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841227) q[0];
sx q[0];
rz(-2.5792482) q[0];
sx q[0];
rz(1.9447295) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9005269) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(-0.17503967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9073346) q[1];
sx q[1];
rz(-1.2163711) q[1];
sx q[1];
rz(-3.1163232) q[1];
rz(0.014724894) q[3];
sx q[3];
rz(-0.71299362) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-3.1299652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25225885) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(2.528119) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83519148) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(-1.4284301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-2.1167397) q[1];
sx q[1];
rz(-3.1139042) q[1];
rz(-2.743268) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(-2.6834965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(-2.7638226) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.4935965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15570116) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(-0.62472384) q[0];
x q[1];
rz(1.1475032) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(2.7763979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6495325) q[1];
sx q[1];
rz(-2.3146555) q[1];
sx q[1];
rz(-0.56708132) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3573523) q[3];
sx q[3];
rz(-2.6297914) q[3];
sx q[3];
rz(-1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.2517694) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(-0.10929508) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(0.072862236) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3423791) q[0];
sx q[0];
rz(-1.8633435) q[0];
sx q[0];
rz(-0.53036687) q[0];
rz(-pi) q[1];
rz(-0.5046919) q[2];
sx q[2];
rz(-2.1345277) q[2];
sx q[2];
rz(1.6187514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2485679) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(-0.47793169) q[1];
rz(-pi) q[2];
rz(-2.7618802) q[3];
sx q[3];
rz(-1.100193) q[3];
sx q[3];
rz(-3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(-0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28829065) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13578829) q[0];
sx q[0];
rz(-3.0514768) q[0];
sx q[0];
rz(-2.3401005) q[0];
rz(0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(-1.136214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.006146487) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(-2.1250493) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0763361) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.4927647) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(-1.5169253) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(-2.0786053) q[2];
sx q[2];
rz(-2.4175736) q[2];
sx q[2];
rz(-2.8833817) q[2];
rz(-1.4528965) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
