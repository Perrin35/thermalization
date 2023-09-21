OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(-0.34531265) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(-2.8023281) q[1];
sx q[1];
rz(-0.27944922) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011443519) q[0];
sx q[0];
rz(-2.7675417) q[0];
sx q[0];
rz(-0.98056294) q[0];
rz(-pi) q[1];
rz(-1.9136393) q[2];
sx q[2];
rz(-2.3228085) q[2];
sx q[2];
rz(2.0057099) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5073519) q[1];
sx q[1];
rz(-1.5841286) q[1];
sx q[1];
rz(1.3807382) q[1];
rz(0.64198288) q[3];
sx q[3];
rz(-1.6025474) q[3];
sx q[3];
rz(0.25105219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7001069) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(2.1842365) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807724) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(1.3445688) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(0.63562524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5168415) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(-1.8952888) q[0];
rz(-pi) q[1];
rz(-1.9102155) q[2];
sx q[2];
rz(-3.0017827) q[2];
sx q[2];
rz(1.8949053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0932514) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(1.3202207) q[1];
x q[2];
rz(-0.73511519) q[3];
sx q[3];
rz(-0.8494091) q[3];
sx q[3];
rz(1.848156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(-1.2403437) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813331) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397944) q[0];
sx q[0];
rz(-2.5520303) q[0];
sx q[0];
rz(-0.88296417) q[0];
rz(-0.065504727) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(0.62578177) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3420521) q[1];
sx q[1];
rz(-0.5040579) q[1];
sx q[1];
rz(-2.6964158) q[1];
rz(-0.86359777) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(1.1532702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(3.0857962) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(-0.24308932) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980804) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(2.9699516) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(2.8780639) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7292273) q[0];
sx q[0];
rz(-2.1365039) q[0];
sx q[0];
rz(1.4501146) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0758923) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(2.2631753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4315223) q[1];
sx q[1];
rz(-2.353851) q[1];
sx q[1];
rz(0.4962994) q[1];
rz(0.27459009) q[3];
sx q[3];
rz(-1.8927186) q[3];
sx q[3];
rz(1.4720751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.8919224) q[2];
rz(-1.194687) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(-0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93976218) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(2.8129541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6304566) q[0];
sx q[0];
rz(-1.6726613) q[0];
sx q[0];
rz(1.0826375) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26942307) q[2];
sx q[2];
rz(-1.7413483) q[2];
sx q[2];
rz(-1.558687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8935834) q[1];
sx q[1];
rz(-1.6343405) q[1];
sx q[1];
rz(0.43339543) q[1];
x q[2];
rz(-2.394649) q[3];
sx q[3];
rz(-1.9455823) q[3];
sx q[3];
rz(2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(0.16061352) q[2];
rz(1.1462071) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(-2.699111) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839358) q[0];
sx q[0];
rz(-1.929495) q[0];
sx q[0];
rz(1.1789807) q[0];
rz(-pi) q[1];
rz(1.2636678) q[2];
sx q[2];
rz(-1.0182292) q[2];
sx q[2];
rz(0.27586684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1861371) q[1];
sx q[1];
rz(-1.326814) q[1];
sx q[1];
rz(-1.2402435) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.093591452) q[3];
sx q[3];
rz(-1.2564661) q[3];
sx q[3];
rz(-2.4623354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.68142146) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(-2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63715315) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-2.1321645) q[0];
rz(0.55074739) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.1632464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79272072) q[0];
sx q[0];
rz(-0.36076818) q[0];
sx q[0];
rz(-1.9804721) q[0];
rz(1.8613569) q[2];
sx q[2];
rz(-0.72939789) q[2];
sx q[2];
rz(1.1682208) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59730676) q[1];
sx q[1];
rz(-1.7339216) q[1];
sx q[1];
rz(0.2225999) q[1];
x q[2];
rz(2.2667415) q[3];
sx q[3];
rz(-1.8112744) q[3];
sx q[3];
rz(1.8868173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(2.3434095) q[2];
rz(-2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(-3.0154862) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(-1.2164446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053596) q[0];
sx q[0];
rz(-1.981712) q[0];
sx q[0];
rz(-2.4634325) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8100796) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(0.31994672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(-1.3575421) q[1];
x q[2];
rz(-1.3835195) q[3];
sx q[3];
rz(-0.98687275) q[3];
sx q[3];
rz(0.54959471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2609666) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(-2.96636) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(-2.9328226) q[0];
rz(-1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-2.0170905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18692423) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(-0.78878553) q[0];
rz(-pi) q[1];
rz(1.8363981) q[2];
sx q[2];
rz(-2.2679066) q[2];
sx q[2];
rz(3.0150989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.323656) q[1];
sx q[1];
rz(-1.7041823) q[1];
sx q[1];
rz(-0.0011841983) q[1];
x q[2];
rz(-2.1379925) q[3];
sx q[3];
rz(-2.3611259) q[3];
sx q[3];
rz(1.6571852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0008529) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(-1.595165) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(-2.3433698) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4231739) q[0];
sx q[0];
rz(-2.383854) q[0];
sx q[0];
rz(2.5916879) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28658861) q[2];
sx q[2];
rz(-1.5550685) q[2];
sx q[2];
rz(-0.77428267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3372911) q[1];
sx q[1];
rz(-1.9807439) q[1];
sx q[1];
rz(-1.0211584) q[1];
rz(2.4622915) q[3];
sx q[3];
rz(-1.7656109) q[3];
sx q[3];
rz(2.4904345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(-0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395441) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(1.7655903) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(-2.3464936) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(1.8134712) q[3];
sx q[3];
rz(-0.96488733) q[3];
sx q[3];
rz(0.78028954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];