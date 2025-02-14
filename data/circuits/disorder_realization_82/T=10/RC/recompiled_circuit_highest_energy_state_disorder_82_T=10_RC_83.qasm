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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(-1.9544344) q[0];
rz(-2.1597593) q[1];
sx q[1];
rz(-2.667006) q[1];
sx q[1];
rz(-2.2003953) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1008074) q[0];
sx q[0];
rz(-1.5691162) q[0];
sx q[0];
rz(0.4953065) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12856483) q[2];
sx q[2];
rz(-1.3701653) q[2];
sx q[2];
rz(-0.085392949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6189673) q[1];
sx q[1];
rz(-1.9272447) q[1];
sx q[1];
rz(0.59774472) q[1];
rz(-3.0676316) q[3];
sx q[3];
rz(-1.1152905) q[3];
sx q[3];
rz(2.2799673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5473951) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(-2.1217864) q[2];
rz(-0.29113537) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(1.3297133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52784598) q[0];
sx q[0];
rz(-2.3044523) q[0];
sx q[0];
rz(-1.6380731) q[0];
rz(1.9095406) q[1];
sx q[1];
rz(-0.82545311) q[1];
sx q[1];
rz(-1.253461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37663662) q[0];
sx q[0];
rz(-0.47237637) q[0];
sx q[0];
rz(-2.0852023) q[0];
rz(2.6297036) q[2];
sx q[2];
rz(-1.8266668) q[2];
sx q[2];
rz(0.40718174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86402363) q[1];
sx q[1];
rz(-1.7146085) q[1];
sx q[1];
rz(0.20419232) q[1];
rz(2.3971629) q[3];
sx q[3];
rz(-1.2334614) q[3];
sx q[3];
rz(-0.51038377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85679179) q[2];
sx q[2];
rz(-0.54232001) q[2];
sx q[2];
rz(-2.8112603) q[2];
rz(3.0387943) q[3];
sx q[3];
rz(-1.2448575) q[3];
sx q[3];
rz(-0.61906329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61515808) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(1.0636299) q[0];
rz(-3.0377153) q[1];
sx q[1];
rz(-0.79236284) q[1];
sx q[1];
rz(2.0361384) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5828667) q[0];
sx q[0];
rz(-1.5387282) q[0];
sx q[0];
rz(-1.523287) q[0];
x q[1];
rz(1.2609405) q[2];
sx q[2];
rz(-0.83854691) q[2];
sx q[2];
rz(-2.1149879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11239219) q[1];
sx q[1];
rz(-1.4489739) q[1];
sx q[1];
rz(1.0531154) q[1];
rz(1.2978992) q[3];
sx q[3];
rz(-3.0981494) q[3];
sx q[3];
rz(-1.3406875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1198472) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(2.9160807) q[2];
rz(1.8146023) q[3];
sx q[3];
rz(-2.3046389) q[3];
sx q[3];
rz(-2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91318146) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(2.5033503) q[0];
rz(2.0300716) q[1];
sx q[1];
rz(-0.94739729) q[1];
sx q[1];
rz(-0.18724719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1094681) q[0];
sx q[0];
rz(-0.95987849) q[0];
sx q[0];
rz(-0.35232589) q[0];
x q[1];
rz(-0.78088607) q[2];
sx q[2];
rz(-1.5263288) q[2];
sx q[2];
rz(-1.5574297) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2350712) q[1];
sx q[1];
rz(-1.491761) q[1];
sx q[1];
rz(2.6071878) q[1];
rz(-pi) q[2];
rz(-1.6583913) q[3];
sx q[3];
rz(-0.85932743) q[3];
sx q[3];
rz(2.8397444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.036666544) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(-0.46910134) q[2];
rz(-0.30096287) q[3];
sx q[3];
rz(-0.260869) q[3];
sx q[3];
rz(1.4891967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.24333532) q[0];
sx q[0];
rz(-1.6247592) q[0];
sx q[0];
rz(0.54018706) q[0];
rz(0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(-0.25925055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40045462) q[0];
sx q[0];
rz(-1.260886) q[0];
sx q[0];
rz(1.8551338) q[0];
rz(-pi) q[1];
rz(-2.0820145) q[2];
sx q[2];
rz(-1.9661742) q[2];
sx q[2];
rz(1.2653654) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.479695) q[1];
sx q[1];
rz(-2.5475218) q[1];
sx q[1];
rz(-1.8516783) q[1];
rz(0.90796538) q[3];
sx q[3];
rz(-2.7304711) q[3];
sx q[3];
rz(-3.079252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79869142) q[2];
sx q[2];
rz(-0.68621245) q[2];
sx q[2];
rz(-0.9101103) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(-2.4250987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27123505) q[0];
sx q[0];
rz(-0.39327455) q[0];
sx q[0];
rz(-1.9140592) q[0];
rz(-0.50499376) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(-1.9379157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6617463) q[0];
sx q[0];
rz(-1.8176314) q[0];
sx q[0];
rz(-1.6936452) q[0];
rz(-2.2199322) q[2];
sx q[2];
rz(-1.1716915) q[2];
sx q[2];
rz(1.6726577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.595049) q[1];
sx q[1];
rz(-2.1814934) q[1];
sx q[1];
rz(-0.86169075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0730482) q[3];
sx q[3];
rz(-2.6274791) q[3];
sx q[3];
rz(2.3069366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(-2.1825979) q[2];
rz(-0.2995019) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(-1.449077) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9129979) q[0];
sx q[0];
rz(-1.8916425) q[0];
sx q[0];
rz(-2.6477497) q[0];
rz(-0.85810581) q[1];
sx q[1];
rz(-1.9789507) q[1];
sx q[1];
rz(-1.7869305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6436967) q[0];
sx q[0];
rz(-0.61031155) q[0];
sx q[0];
rz(-1.1417139) q[0];
rz(1.0787215) q[2];
sx q[2];
rz(-1.4495696) q[2];
sx q[2];
rz(1.4722784) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2676182) q[1];
sx q[1];
rz(-0.81470352) q[1];
sx q[1];
rz(1.6674629) q[1];
rz(2.1936761) q[3];
sx q[3];
rz(-2.3661032) q[3];
sx q[3];
rz(2.754107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6845503) q[2];
sx q[2];
rz(-0.6251308) q[2];
sx q[2];
rz(-0.54614145) q[2];
rz(-0.1977194) q[3];
sx q[3];
rz(-2.1106014) q[3];
sx q[3];
rz(2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2208743) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(-2.0350631) q[0];
rz(0.58440343) q[1];
sx q[1];
rz(-1.1738651) q[1];
sx q[1];
rz(2.1388785) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4179392) q[0];
sx q[0];
rz(-2.0040211) q[0];
sx q[0];
rz(2.9648925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2539576) q[2];
sx q[2];
rz(-2.7681367) q[2];
sx q[2];
rz(1.4591726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3467873) q[1];
sx q[1];
rz(-2.0827939) q[1];
sx q[1];
rz(-1.6122628) q[1];
rz(-pi) q[2];
rz(-1.2190621) q[3];
sx q[3];
rz(-2.1603843) q[3];
sx q[3];
rz(1.0603348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5661085) q[2];
sx q[2];
rz(-2.2169952) q[2];
sx q[2];
rz(1.7522579) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(-2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.093512) q[0];
sx q[0];
rz(-2.7219613) q[0];
sx q[0];
rz(-2.804948) q[0];
rz(3.0632784) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(-1.1557109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0037831941) q[0];
sx q[0];
rz(-2.3562618) q[0];
sx q[0];
rz(2.0565874) q[0];
x q[1];
rz(-1.1572273) q[2];
sx q[2];
rz(-1.509827) q[2];
sx q[2];
rz(0.79198906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48574572) q[1];
sx q[1];
rz(-2.0585723) q[1];
sx q[1];
rz(-2.9887245) q[1];
x q[2];
rz(-1.5404245) q[3];
sx q[3];
rz(-1.8874468) q[3];
sx q[3];
rz(0.79398549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.081850514) q[2];
sx q[2];
rz(-2.0589224) q[2];
sx q[2];
rz(2.3183863) q[2];
rz(-2.1246223) q[3];
sx q[3];
rz(-2.7401994) q[3];
sx q[3];
rz(-0.038173525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80242771) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(-1.4270576) q[0];
rz(1.4786221) q[1];
sx q[1];
rz(-1.0717816) q[1];
sx q[1];
rz(0.85085416) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6991147) q[0];
sx q[0];
rz(-0.48865396) q[0];
sx q[0];
rz(1.3241934) q[0];
rz(-1.698174) q[2];
sx q[2];
rz(-1.4847404) q[2];
sx q[2];
rz(2.0620777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58356279) q[1];
sx q[1];
rz(-1.0475912) q[1];
sx q[1];
rz(-1.5091404) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6353278) q[3];
sx q[3];
rz(-1.7108265) q[3];
sx q[3];
rz(-2.9794995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.458638) q[2];
sx q[2];
rz(-1.3502324) q[2];
sx q[2];
rz(-0.82009912) q[2];
rz(1.594515) q[3];
sx q[3];
rz(-1.4328512) q[3];
sx q[3];
rz(1.9471751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1632669) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(-2.3691879) q[1];
sx q[1];
rz(-0.62320566) q[1];
sx q[1];
rz(-2.1355245) q[1];
rz(-1.7987235) q[2];
sx q[2];
rz(-0.89222903) q[2];
sx q[2];
rz(-1.0747838) q[2];
rz(2.0396654) q[3];
sx q[3];
rz(-2.0004604) q[3];
sx q[3];
rz(-0.018044005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
