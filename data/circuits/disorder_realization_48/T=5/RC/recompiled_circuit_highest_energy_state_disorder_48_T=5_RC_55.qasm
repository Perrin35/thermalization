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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(-1.7755427) q[1];
sx q[1];
rz(-0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1695367) q[0];
sx q[0];
rz(-2.0436624) q[0];
sx q[0];
rz(-2.1465143) q[0];
x q[1];
rz(1.0452966) q[2];
sx q[2];
rz(-2.3901148) q[2];
sx q[2];
rz(0.95096248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9626689) q[1];
sx q[1];
rz(-1.0826063) q[1];
sx q[1];
rz(2.8250384) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3955807) q[3];
sx q[3];
rz(-1.0300582) q[3];
sx q[3];
rz(-2.9313425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-1.0203699) q[2];
rz(-0.72888914) q[3];
sx q[3];
rz(-2.708669) q[3];
sx q[3];
rz(-1.5322022) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6272524) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(0.038272055) q[0];
rz(-0.12380869) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(-1.5708057) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8077652) q[0];
sx q[0];
rz(-3.1196731) q[0];
sx q[0];
rz(-2.0448565) q[0];
x q[1];
rz(0.77124243) q[2];
sx q[2];
rz(-1.6420157) q[2];
sx q[2];
rz(2.3598755) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99421895) q[1];
sx q[1];
rz(-0.0011073907) q[1];
sx q[1];
rz(1.1785517) q[1];
rz(-pi) q[2];
rz(2.6255076) q[3];
sx q[3];
rz(-1.4654963) q[3];
sx q[3];
rz(-0.95296958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(0.5788571) q[2];
rz(-1.045643) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(0.75700179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-2.6876167) q[0];
rz(0.16481915) q[1];
sx q[1];
rz(-0.77607981) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904083) q[0];
sx q[0];
rz(-1.7421573) q[0];
sx q[0];
rz(1.1771855) q[0];
rz(-pi) q[1];
rz(1.0452095) q[2];
sx q[2];
rz(-1.8586577) q[2];
sx q[2];
rz(2.9381816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3537703) q[1];
sx q[1];
rz(-0.90836891) q[1];
sx q[1];
rz(1.1819928) q[1];
rz(-pi) q[2];
rz(-1.0461058) q[3];
sx q[3];
rz(-0.98011049) q[3];
sx q[3];
rz(-2.5928796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(-1.3564159) q[2];
rz(0.49946579) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(-2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.22719638) q[0];
sx q[0];
rz(-2.2375186) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.1610441) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6446021) q[0];
sx q[0];
rz(-1.0287026) q[0];
sx q[0];
rz(-1.5250456) q[0];
x q[1];
rz(-0.19660825) q[2];
sx q[2];
rz(-0.93485445) q[2];
sx q[2];
rz(-1.3351118) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30018729) q[1];
sx q[1];
rz(-2.3267728) q[1];
sx q[1];
rz(0.38736613) q[1];
rz(-pi) q[2];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78493541) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.7735205) q[2];
rz(1.8562227) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042498978) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(-1.7294783) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(1.5589421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928685) q[0];
sx q[0];
rz(-0.72843116) q[0];
sx q[0];
rz(0.65653519) q[0];
rz(-pi) q[1];
rz(-0.71835235) q[2];
sx q[2];
rz(-1.5081154) q[2];
sx q[2];
rz(2.2557392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0274899) q[1];
sx q[1];
rz(-1.5943502) q[1];
sx q[1];
rz(-1.0204169) q[1];
rz(1.8924501) q[3];
sx q[3];
rz(-1.6603004) q[3];
sx q[3];
rz(-2.9200875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5530508) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2734964) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-2.9582276) q[0];
rz(-0.49194899) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(1.2194182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6693319) q[0];
sx q[0];
rz(-0.12024256) q[0];
sx q[0];
rz(0.73445102) q[0];
rz(-1.2017131) q[2];
sx q[2];
rz(-0.08412349) q[2];
sx q[2];
rz(-0.28757986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9907836) q[1];
sx q[1];
rz(-1.1649302) q[1];
sx q[1];
rz(1.1836786) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0571396) q[3];
sx q[3];
rz(-1.1005529) q[3];
sx q[3];
rz(-0.081428226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(-2.2149337) q[3];
sx q[3];
rz(-2.4121273) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6019186) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(-2.4251921) q[0];
rz(2.7377103) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(-1.7049888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88750171) q[0];
sx q[0];
rz(-1.1932826) q[0];
sx q[0];
rz(0.15477263) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.613343) q[2];
sx q[2];
rz(-1.4563603) q[2];
sx q[2];
rz(-0.50811646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75175367) q[1];
sx q[1];
rz(-1.1734278) q[1];
sx q[1];
rz(1.6786871) q[1];
rz(-2.3327504) q[3];
sx q[3];
rz(-1.0029965) q[3];
sx q[3];
rz(2.0140935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7776103) q[2];
sx q[2];
rz(-0.97475514) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(-0.86723793) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(-1.3245827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6548178) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(0.73100334) q[0];
rz(-0.77516088) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(-2.7269272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054798445) q[0];
sx q[0];
rz(-2.2761184) q[0];
sx q[0];
rz(0.84673015) q[0];
x q[1];
rz(-2.6111905) q[2];
sx q[2];
rz(-1.7549522) q[2];
sx q[2];
rz(3.0006486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4463628) q[1];
sx q[1];
rz(-1.8703504) q[1];
sx q[1];
rz(1.8336373) q[1];
x q[2];
rz(-3.0075668) q[3];
sx q[3];
rz(-1.408159) q[3];
sx q[3];
rz(-2.7664879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23184648) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(-0.20364729) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(-0.017729433) q[0];
rz(2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(-2.0808992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6700376) q[0];
sx q[0];
rz(-3.0714543) q[0];
sx q[0];
rz(0.20151968) q[0];
x q[1];
rz(-2.1589375) q[2];
sx q[2];
rz(-1.3545631) q[2];
sx q[2];
rz(-0.71046358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0073787) q[1];
sx q[1];
rz(-1.7430647) q[1];
sx q[1];
rz(-0.80165792) q[1];
x q[2];
rz(0.10194998) q[3];
sx q[3];
rz(-0.52604874) q[3];
sx q[3];
rz(-2.7395484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1549418) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(0.6443392) q[2];
rz(-0.84195697) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7613206) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(-0.45183387) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(-1.5712646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1399978) q[0];
sx q[0];
rz(-0.17029914) q[0];
sx q[0];
rz(0.036080555) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0250835) q[2];
sx q[2];
rz(-1.0986058) q[2];
sx q[2];
rz(1.0509059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3305407) q[1];
sx q[1];
rz(-1.1833541) q[1];
sx q[1];
rz(0.13918332) q[1];
rz(-pi) q[2];
rz(0.7223981) q[3];
sx q[3];
rz(-1.5766451) q[3];
sx q[3];
rz(1.7286144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(-0.10562854) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-2.6251729) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.327772) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(-2.9227921) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(-1.1458746) q[2];
sx q[2];
rz(-2.6617202) q[2];
sx q[2];
rz(-0.94750994) q[2];
rz(0.65746751) q[3];
sx q[3];
rz(-1.8651265) q[3];
sx q[3];
rz(-0.65262564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
