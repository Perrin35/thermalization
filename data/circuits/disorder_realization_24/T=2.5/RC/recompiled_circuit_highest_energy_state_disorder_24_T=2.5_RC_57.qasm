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
rz(2.6296122) q[0];
sx q[0];
rz(-0.57096243) q[0];
sx q[0];
rz(2.9538739) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(-1.6698281) q[1];
sx q[1];
rz(-1.8522813) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780802) q[0];
sx q[0];
rz(-0.31766444) q[0];
sx q[0];
rz(-2.4207522) q[0];
rz(0.046229049) q[2];
sx q[2];
rz(-2.5704489) q[2];
sx q[2];
rz(-2.2276304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.047886176) q[1];
sx q[1];
rz(-2.0290861) q[1];
sx q[1];
rz(1.8480193) q[1];
x q[2];
rz(-2.1407152) q[3];
sx q[3];
rz(-0.412768) q[3];
sx q[3];
rz(-1.643553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(2.3647986) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(-1.2031901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5341107) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(2.7718995) q[1];
sx q[1];
rz(-2.1001215) q[1];
sx q[1];
rz(2.0236156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48644439) q[0];
sx q[0];
rz(-2.110743) q[0];
sx q[0];
rz(-2.6686431) q[0];
rz(-pi) q[1];
rz(2.4547365) q[2];
sx q[2];
rz(-2.0389028) q[2];
sx q[2];
rz(1.3324225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4395968) q[1];
sx q[1];
rz(-1.1267822) q[1];
sx q[1];
rz(0.16298144) q[1];
rz(-pi) q[2];
rz(0.93922521) q[3];
sx q[3];
rz(-1.3520006) q[3];
sx q[3];
rz(-0.18635264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2758241) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(-2.7562874) q[2];
rz(0.038711874) q[3];
sx q[3];
rz(-0.35274115) q[3];
sx q[3];
rz(-0.64046162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63967079) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(-1.84024) q[0];
rz(-3.0838857) q[1];
sx q[1];
rz(-2.6951908) q[1];
sx q[1];
rz(1.8147963) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037687) q[0];
sx q[0];
rz(-0.7531727) q[0];
sx q[0];
rz(-2.2979548) q[0];
rz(-pi) q[1];
rz(-2.9811601) q[2];
sx q[2];
rz(-1.7084873) q[2];
sx q[2];
rz(1.4169803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51242764) q[1];
sx q[1];
rz(-1.8574134) q[1];
sx q[1];
rz(-2.5932339) q[1];
x q[2];
rz(1.9671828) q[3];
sx q[3];
rz(-0.76226888) q[3];
sx q[3];
rz(-1.1571194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3441299) q[2];
sx q[2];
rz(-2.3064488) q[2];
sx q[2];
rz(-0.75378913) q[2];
rz(-1.3695184) q[3];
sx q[3];
rz(-1.4621719) q[3];
sx q[3];
rz(2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(-1.1766379) q[1];
sx q[1];
rz(-1.6329012) q[1];
sx q[1];
rz(-0.36697695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4876539) q[0];
sx q[0];
rz(-1.5040845) q[0];
sx q[0];
rz(-0.056746527) q[0];
x q[1];
rz(2.4796333) q[2];
sx q[2];
rz(-2.7267704) q[2];
sx q[2];
rz(-2.7560459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5287197) q[1];
sx q[1];
rz(-2.7320128) q[1];
sx q[1];
rz(0.50480857) q[1];
rz(-2.9212679) q[3];
sx q[3];
rz(-0.37543618) q[3];
sx q[3];
rz(2.8130949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(1.270594) q[2];
rz(-2.1805084) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(-3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56080317) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(-1.8286937) q[0];
rz(-1.987223) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(0.55783522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9168222) q[0];
sx q[0];
rz(-1.0616515) q[0];
sx q[0];
rz(-0.96352412) q[0];
rz(-pi) q[1];
rz(1.6480321) q[2];
sx q[2];
rz(-1.4106515) q[2];
sx q[2];
rz(0.84116615) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.085943) q[1];
sx q[1];
rz(-2.331877) q[1];
sx q[1];
rz(1.9941866) q[1];
rz(-pi) q[2];
rz(1.722103) q[3];
sx q[3];
rz(-2.0278721) q[3];
sx q[3];
rz(-0.099927038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(0.84929973) q[2];
rz(-0.15527209) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.3575386) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(-2.6724755) q[0];
rz(-0.47438374) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(-2.3862086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635272) q[0];
sx q[0];
rz(-1.2995412) q[0];
sx q[0];
rz(3.1370107) q[0];
x q[1];
rz(-2.7003661) q[2];
sx q[2];
rz(-0.6775113) q[2];
sx q[2];
rz(2.4063619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1052221) q[1];
sx q[1];
rz(-1.3579988) q[1];
sx q[1];
rz(-0.12963055) q[1];
rz(-pi) q[2];
rz(2.3085318) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(1.2556094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0581806) q[2];
sx q[2];
rz(-1.8184793) q[2];
sx q[2];
rz(-0.18784909) q[2];
rz(1.2457054) q[3];
sx q[3];
rz(-1.7626423) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.93194) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(2.7506822) q[0];
rz(0.93005013) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(1.3040868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025607312) q[0];
sx q[0];
rz(-1.5497594) q[0];
sx q[0];
rz(1.2198835) q[0];
rz(-pi) q[1];
rz(-1.6565336) q[2];
sx q[2];
rz(-1.9112327) q[2];
sx q[2];
rz(2.1742333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2432989) q[1];
sx q[1];
rz(-1.1203246) q[1];
sx q[1];
rz(-0.18646481) q[1];
rz(-pi) q[2];
rz(1.5518858) q[3];
sx q[3];
rz(-1.2745556) q[3];
sx q[3];
rz(-0.66413524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2685252) q[2];
sx q[2];
rz(-0.26669845) q[2];
sx q[2];
rz(-0.11338691) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(-2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78405821) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(-0.12406021) q[0];
rz(0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(-1.6990936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47114357) q[0];
sx q[0];
rz(-1.6713977) q[0];
sx q[0];
rz(1.1858245) q[0];
rz(0.14367468) q[2];
sx q[2];
rz(-0.29479237) q[2];
sx q[2];
rz(2.2058892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6119163) q[1];
sx q[1];
rz(-1.2280591) q[1];
sx q[1];
rz(-0.14118282) q[1];
rz(-pi) q[2];
rz(0.7300709) q[3];
sx q[3];
rz(-1.7582809) q[3];
sx q[3];
rz(-1.2227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1759935) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(-2.3967801) q[2];
rz(-2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-2.9162245) q[0];
rz(1.1124181) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(-0.89881277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3254125) q[0];
sx q[0];
rz(-2.1006148) q[0];
sx q[0];
rz(0.90318944) q[0];
x q[1];
rz(0.5139047) q[2];
sx q[2];
rz(-1.9348839) q[2];
sx q[2];
rz(-0.68825005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1370217) q[1];
sx q[1];
rz(-0.70024509) q[1];
sx q[1];
rz(1.2546468) q[1];
rz(-pi) q[2];
rz(-3.131116) q[3];
sx q[3];
rz(-1.8636522) q[3];
sx q[3];
rz(-1.7550857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-2.782235) q[2];
rz(-0.25092956) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(-0.8771483) q[0];
rz(1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(2.3559779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025185275) q[0];
sx q[0];
rz(-1.3337564) q[0];
sx q[0];
rz(0.22444207) q[0];
x q[1];
rz(2.4841275) q[2];
sx q[2];
rz(-0.6415002) q[2];
sx q[2];
rz(-0.408022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12745107) q[1];
sx q[1];
rz(-1.9147062) q[1];
sx q[1];
rz(-1.6456855) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26667554) q[3];
sx q[3];
rz(-0.52191041) q[3];
sx q[3];
rz(0.44346186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(-0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(2.8452528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7947163) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(0.90732668) q[1];
sx q[1];
rz(-1.0846039) q[1];
sx q[1];
rz(-0.32122282) q[1];
rz(-0.39225929) q[2];
sx q[2];
rz(-1.0259368) q[2];
sx q[2];
rz(-2.1403014) q[2];
rz(2.6014544) q[3];
sx q[3];
rz(-2.028699) q[3];
sx q[3];
rz(1.5069458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
