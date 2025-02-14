OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.59392053) q[0];
sx q[0];
rz(-2.5424818) q[0];
sx q[0];
rz(0.17679086) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(3.6511753) q[1];
sx q[1];
rz(9.3378172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158543) q[0];
sx q[0];
rz(-0.70924067) q[0];
sx q[0];
rz(-2.5410557) q[0];
x q[1];
rz(0.571351) q[2];
sx q[2];
rz(-1.4077912) q[2];
sx q[2];
rz(-0.75955078) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46087825) q[1];
sx q[1];
rz(-1.6986004) q[1];
sx q[1];
rz(3.0640814) q[1];
x q[2];
rz(-2.7417408) q[3];
sx q[3];
rz(-0.65910554) q[3];
sx q[3];
rz(-0.53149283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1842492) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(2.2653968) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0210719) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(0.72129321) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(0.24478197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1230303) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(0.91482343) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3201687) q[2];
sx q[2];
rz(-1.9000179) q[2];
sx q[2];
rz(-0.034857817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2657425) q[1];
sx q[1];
rz(-2.2835554) q[1];
sx q[1];
rz(-0.73709388) q[1];
x q[2];
rz(1.7177204) q[3];
sx q[3];
rz(-1.8645446) q[3];
sx q[3];
rz(2.5054641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038079809) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(0.14347759) q[0];
rz(1.4682651) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(-3.0751244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398427) q[0];
sx q[0];
rz(-0.89213138) q[0];
sx q[0];
rz(-1.6854273) q[0];
x q[1];
rz(1.1622381) q[2];
sx q[2];
rz(-2.5600932) q[2];
sx q[2];
rz(1.3571908) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4926949) q[1];
sx q[1];
rz(-1.4543415) q[1];
sx q[1];
rz(-1.4728738) q[1];
x q[2];
rz(-1.4044365) q[3];
sx q[3];
rz(-2.3041953) q[3];
sx q[3];
rz(-0.93033965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1116144) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(2.8800268) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12544352) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(0.087015986) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(-3.0932025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3589879) q[0];
sx q[0];
rz(-1.6986956) q[0];
sx q[0];
rz(2.3652698) q[0];
rz(-2.9913285) q[2];
sx q[2];
rz(-1.3925526) q[2];
sx q[2];
rz(1.5706289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.794487) q[1];
sx q[1];
rz(-2.2655267) q[1];
sx q[1];
rz(-1.6415651) q[1];
rz(-pi) q[2];
rz(-2.5278306) q[3];
sx q[3];
rz(-2.1832972) q[3];
sx q[3];
rz(2.4341754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(0.44130138) q[2];
rz(-0.86853164) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(1.3335479) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(1.6983039) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.7194933) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.559645) q[0];
sx q[0];
rz(-1.0302246) q[0];
sx q[0];
rz(-1.7339043) q[0];
rz(-0.70962064) q[2];
sx q[2];
rz(-0.31236744) q[2];
sx q[2];
rz(-1.4210977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15502587) q[1];
sx q[1];
rz(-1.4490713) q[1];
sx q[1];
rz(-1.6568068) q[1];
rz(-pi) q[2];
rz(0.57797076) q[3];
sx q[3];
rz(-2.1506566) q[3];
sx q[3];
rz(2.487929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(1.5606073) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(-0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.018589858) q[0];
sx q[0];
rz(-0.81273166) q[0];
sx q[0];
rz(2.4216968) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(0.24615157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6841078) q[0];
sx q[0];
rz(-0.2817758) q[0];
sx q[0];
rz(-3.0578567) q[0];
x q[1];
rz(-2.7163137) q[2];
sx q[2];
rz(-2.3375247) q[2];
sx q[2];
rz(-0.67415392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64061058) q[1];
sx q[1];
rz(-0.74509186) q[1];
sx q[1];
rz(0.19265811) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6443374) q[3];
sx q[3];
rz(-0.75102931) q[3];
sx q[3];
rz(-2.8145144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2972357) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(-0.921184) q[0];
rz(1.4211897) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9807905) q[0];
sx q[0];
rz(-1.9483999) q[0];
sx q[0];
rz(2.9016205) q[0];
x q[1];
rz(-0.065804577) q[2];
sx q[2];
rz(-1.9901681) q[2];
sx q[2];
rz(2.9247627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91749633) q[1];
sx q[1];
rz(-2.4451548) q[1];
sx q[1];
rz(0.68989474) q[1];
x q[2];
rz(2.7217676) q[3];
sx q[3];
rz(-1.9732765) q[3];
sx q[3];
rz(-0.67129204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(2.7040238) q[2];
rz(0.82018745) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95847982) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(2.130765) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6626604) q[0];
sx q[0];
rz(-0.17864431) q[0];
sx q[0];
rz(-1.7500072) q[0];
x q[1];
rz(-2.5055614) q[2];
sx q[2];
rz(-2.6205728) q[2];
sx q[2];
rz(2.4907559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6514142) q[1];
sx q[1];
rz(-1.7140577) q[1];
sx q[1];
rz(2.4552022) q[1];
rz(-0.88415481) q[3];
sx q[3];
rz(-3.0100689) q[3];
sx q[3];
rz(-2.2145074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7944472) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(-2.6084206) q[2];
rz(1.0779856) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(2.3593498) q[0];
rz(3.0714463) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(2.9152962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057581456) q[0];
sx q[0];
rz(-1.6414483) q[0];
sx q[0];
rz(-1.6397301) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3391206) q[2];
sx q[2];
rz(-2.4501928) q[2];
sx q[2];
rz(-2.0350128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28412425) q[1];
sx q[1];
rz(-2.5924666) q[1];
sx q[1];
rz(-0.21128564) q[1];
rz(-pi) q[2];
rz(-1.3198648) q[3];
sx q[3];
rz(-1.5490412) q[3];
sx q[3];
rz(1.3298498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1066863) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(0.26257026) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(-3.0098651) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866289) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(-0.18375272) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(0.62350887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348891) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(1.2280653) q[0];
rz(-0.76304014) q[2];
sx q[2];
rz(-2.3614778) q[2];
sx q[2];
rz(0.060465079) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7694313) q[1];
sx q[1];
rz(-1.6767394) q[1];
sx q[1];
rz(1.547936) q[1];
rz(-1.4922114) q[3];
sx q[3];
rz(-2.7659263) q[3];
sx q[3];
rz(-0.54822719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(3.1039216) q[2];
rz(-0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2035718) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(2.6899295) q[1];
sx q[1];
rz(-1.3126806) q[1];
sx q[1];
rz(-1.5246593) q[1];
rz(1.8658609) q[2];
sx q[2];
rz(-1.0574592) q[2];
sx q[2];
rz(-1.9707373) q[2];
rz(0.59521159) q[3];
sx q[3];
rz(-0.32929904) q[3];
sx q[3];
rz(2.3198447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
