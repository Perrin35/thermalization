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
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(-0.018215608) q[0];
rz(-pi) q[1];
rz(-0.78975241) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(1.5577158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5684811) q[1];
sx q[1];
rz(-0.1460349) q[1];
sx q[1];
rz(-2.5558429) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53332897) q[3];
sx q[3];
rz(-1.1733857) q[3];
sx q[3];
rz(1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-2.1437936) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.6404023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.507507) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(0.12269845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31200774) q[2];
sx q[2];
rz(-2.0199676) q[2];
sx q[2];
rz(-2.1829407) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7306819) q[1];
sx q[1];
rz(-2.4446179) q[1];
sx q[1];
rz(2.8721786) q[1];
rz(1.5562015) q[3];
sx q[3];
rz(-2.7505034) q[3];
sx q[3];
rz(-0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.6378145) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2730946) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(-0.065386535) q[0];
x q[1];
rz(2.1554699) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(-3.0509146) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7392002) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(1.6154358) q[1];
rz(-pi) q[2];
rz(1.4012666) q[3];
sx q[3];
rz(-2.4781514) q[3];
sx q[3];
rz(1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(0.7437931) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(-1.9925041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.8115037) q[0];
sx q[0];
rz(1.4562777) q[0];
x q[1];
rz(-0.43196584) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(-0.84749046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60050868) q[1];
sx q[1];
rz(-0.87227548) q[1];
sx q[1];
rz(2.7545616) q[1];
rz(0.37255128) q[3];
sx q[3];
rz(-0.77790341) q[3];
sx q[3];
rz(-2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887061) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-1.0523798) q[0];
rz(-1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(2.2573684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0852172) q[0];
sx q[0];
rz(-1.3803122) q[0];
sx q[0];
rz(-2.9604244) q[0];
rz(-2.7344633) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(-1.7083573) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21806949) q[1];
sx q[1];
rz(-1.4265718) q[1];
sx q[1];
rz(1.5037392) q[1];
x q[2];
rz(2.1518425) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(2.1637672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-0.75105914) q[3];
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
rz(1.7364175) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-2.172519) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5613865) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(1.4617306) q[0];
rz(1.5029807) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(-0.80578795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5507257) q[1];
sx q[1];
rz(-2.496965) q[1];
sx q[1];
rz(2.4025737) q[1];
rz(2.6763775) q[3];
sx q[3];
rz(-2.0485268) q[3];
sx q[3];
rz(-0.623869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(-2.7039841) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.70599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3754417) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(0.96813162) q[0];
rz(0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(-2.6798623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0457397) q[1];
sx q[1];
rz(-2.6624655) q[1];
sx q[1];
rz(2.2687001) q[1];
rz(0.42922677) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(2.6847117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5196461) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(-1.7072862) q[0];
rz(-pi) q[1];
rz(-1.6985967) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(0.70805659) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9284968) q[1];
sx q[1];
rz(-2.1017385) q[1];
sx q[1];
rz(-1.5099768) q[1];
rz(0.41994628) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7565219) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(-0.71682741) q[0];
rz(-pi) q[1];
rz(-2.6529979) q[2];
sx q[2];
rz(-1.5602419) q[2];
sx q[2];
rz(-0.34305629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98665806) q[1];
sx q[1];
rz(-2.0866052) q[1];
sx q[1];
rz(0.45706473) q[1];
rz(-2.4273848) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
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
rz(-3.0420619) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-0.038169233) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(1.5241745) q[0];
rz(-pi) q[1];
rz(1.3456887) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(0.12347808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1306418) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.6846659) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60572259) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(3.0063418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(2.9057124) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.11001982) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(-2.2723972) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
