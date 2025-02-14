OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(7.5543348) q[0];
sx q[0];
rz(13.332097) q[0];
rz(-0.39628059) q[1];
sx q[1];
rz(-0.09274617) q[1];
sx q[1];
rz(0.5527817) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39550135) q[0];
sx q[0];
rz(-1.9736909) q[0];
sx q[0];
rz(2.7159116) q[0];
rz(0.67019318) q[2];
sx q[2];
rz(-2.2171221) q[2];
sx q[2];
rz(1.0431511) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1172792) q[1];
sx q[1];
rz(-1.8872042) q[1];
sx q[1];
rz(-0.98543075) q[1];
rz(-pi) q[2];
rz(1.0200653) q[3];
sx q[3];
rz(-2.4006776) q[3];
sx q[3];
rz(1.7307841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7736241) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(-3.0527414) q[2];
rz(0.39696524) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(-1.4575492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001584) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(2.2870824) q[0];
rz(-0.62141934) q[1];
sx q[1];
rz(-2.8050551) q[1];
sx q[1];
rz(-1.052676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5993505) q[0];
sx q[0];
rz(-2.9711038) q[0];
sx q[0];
rz(-0.64298274) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82727716) q[2];
sx q[2];
rz(-0.57663871) q[2];
sx q[2];
rz(2.1183543) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4667644) q[1];
sx q[1];
rz(-2.383814) q[1];
sx q[1];
rz(1.8653052) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0016453513) q[3];
sx q[3];
rz(-2.1747394) q[3];
sx q[3];
rz(1.920948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66775995) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(0.94998002) q[2];
rz(2.1330323) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(-0.27957255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94473332) q[0];
sx q[0];
rz(-1.2161398) q[0];
sx q[0];
rz(1.3470294) q[0];
rz(-1.0579717) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(3.0314441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63907097) q[0];
sx q[0];
rz(-0.57146954) q[0];
sx q[0];
rz(2.1854464) q[0];
rz(-pi) q[1];
rz(1.7188495) q[2];
sx q[2];
rz(-1.4880991) q[2];
sx q[2];
rz(-0.18872866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85903214) q[1];
sx q[1];
rz(-0.60884005) q[1];
sx q[1];
rz(0.18688272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0581117) q[3];
sx q[3];
rz(-1.1705453) q[3];
sx q[3];
rz(2.6292173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18016732) q[2];
sx q[2];
rz(-1.5632997) q[2];
sx q[2];
rz(-0.041672826) q[2];
rz(2.9336119) q[3];
sx q[3];
rz(-2.9198923) q[3];
sx q[3];
rz(-0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792849) q[0];
sx q[0];
rz(-1.3212181) q[0];
sx q[0];
rz(0.98096171) q[0];
rz(-0.43308577) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(1.513419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5034803) q[0];
sx q[0];
rz(-2.7013218) q[0];
sx q[0];
rz(-2.7087337) q[0];
rz(-pi) q[1];
rz(0.01362205) q[2];
sx q[2];
rz(-0.8994461) q[2];
sx q[2];
rz(-1.7126132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1867552) q[1];
sx q[1];
rz(-1.673497) q[1];
sx q[1];
rz(-2.0689194) q[1];
rz(-pi) q[2];
rz(-1.1918729) q[3];
sx q[3];
rz(-0.66693587) q[3];
sx q[3];
rz(2.2569424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5341586) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(0.81238166) q[2];
rz(1.8937998) q[3];
sx q[3];
rz(-1.2500074) q[3];
sx q[3];
rz(-1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-1.0220818) q[0];
sx q[0];
rz(1.6988276) q[0];
rz(-0.45817786) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(0.75622574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6578015) q[0];
sx q[0];
rz(-1.4157802) q[0];
sx q[0];
rz(2.20138) q[0];
x q[1];
rz(2.5014072) q[2];
sx q[2];
rz(-1.1220896) q[2];
sx q[2];
rz(-1.9404836) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9994558) q[1];
sx q[1];
rz(-0.161006) q[1];
sx q[1];
rz(-0.57438897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4706572) q[3];
sx q[3];
rz(-2.7960145) q[3];
sx q[3];
rz(2.341193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68236399) q[2];
sx q[2];
rz(-1.7533147) q[2];
sx q[2];
rz(2.1759822) q[2];
rz(-0.56973488) q[3];
sx q[3];
rz(-1.1252879) q[3];
sx q[3];
rz(1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4023034) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(1.385561) q[0];
rz(0.7631453) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(3.0398583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67416009) q[0];
sx q[0];
rz(-1.4786353) q[0];
sx q[0];
rz(-2.1919495) q[0];
rz(3.0252803) q[2];
sx q[2];
rz(-0.38517932) q[2];
sx q[2];
rz(2.0907837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6468874) q[1];
sx q[1];
rz(-0.27783074) q[1];
sx q[1];
rz(1.7956514) q[1];
x q[2];
rz(-2.2564628) q[3];
sx q[3];
rz(-0.99381522) q[3];
sx q[3];
rz(2.5654405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0330641) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(2.4326883) q[2];
rz(-2.3809643) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(-0.93543783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.075542299) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(0.72108889) q[0];
rz(-0.9779633) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.579938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1051416) q[0];
sx q[0];
rz(-1.4976964) q[0];
sx q[0];
rz(-0.076687254) q[0];
x q[1];
rz(0.97550895) q[2];
sx q[2];
rz(-1.7726608) q[2];
sx q[2];
rz(-1.6109811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5389293) q[1];
sx q[1];
rz(-2.437535) q[1];
sx q[1];
rz(0.62967794) q[1];
rz(2.0098814) q[3];
sx q[3];
rz(-0.71936047) q[3];
sx q[3];
rz(2.2124825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(1.52012) q[2];
rz(0.4367477) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(-2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590416) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(2.3178597) q[0];
rz(1.1388904) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(1.1762071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0498031) q[0];
sx q[0];
rz(-0.50423008) q[0];
sx q[0];
rz(-1.4780028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61061956) q[2];
sx q[2];
rz(-1.23252) q[2];
sx q[2];
rz(-0.18649292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8835134) q[1];
sx q[1];
rz(-0.23830676) q[1];
sx q[1];
rz(-0.041447354) q[1];
rz(1.6105936) q[3];
sx q[3];
rz(-2.1615513) q[3];
sx q[3];
rz(-1.9692957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1463683) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(2.5541019) q[2];
rz(-0.38582173) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189482) q[0];
sx q[0];
rz(-2.1653439) q[0];
sx q[0];
rz(-2.5318085) q[0];
rz(1.7976286) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(2.9885898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55930821) q[0];
sx q[0];
rz(-2.0366208) q[0];
sx q[0];
rz(1.5488173) q[0];
x q[1];
rz(1.581089) q[2];
sx q[2];
rz(-2.3659945) q[2];
sx q[2];
rz(0.44558172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3045522) q[1];
sx q[1];
rz(-2.4992832) q[1];
sx q[1];
rz(1.539597) q[1];
x q[2];
rz(0.47346799) q[3];
sx q[3];
rz(-0.93515474) q[3];
sx q[3];
rz(1.9994232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1028221) q[2];
sx q[2];
rz(-1.5615347) q[2];
sx q[2];
rz(-2.5193396) q[2];
rz(1.2004987) q[3];
sx q[3];
rz(-1.7222722) q[3];
sx q[3];
rz(-0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(2.0390873) q[1];
sx q[1];
rz(-0.016977221) q[1];
sx q[1];
rz(2.4818518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44910494) q[0];
sx q[0];
rz(-0.67506719) q[0];
sx q[0];
rz(-2.8646742) q[0];
x q[1];
rz(-3.0142886) q[2];
sx q[2];
rz(-1.1603171) q[2];
sx q[2];
rz(2.8212027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28478947) q[1];
sx q[1];
rz(-1.2592013) q[1];
sx q[1];
rz(1.2272528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2485178) q[3];
sx q[3];
rz(-1.5297946) q[3];
sx q[3];
rz(1.463691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3416662) q[2];
sx q[2];
rz(-3.0604) q[2];
sx q[2];
rz(-2.9426835) q[2];
rz(-0.79814664) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(1.859349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2220919) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(0.78492289) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(1.2047819) q[2];
sx q[2];
rz(-0.13199619) q[2];
sx q[2];
rz(0.97313626) q[2];
rz(-1.296923) q[3];
sx q[3];
rz(-0.92007888) q[3];
sx q[3];
rz(-2.32213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
