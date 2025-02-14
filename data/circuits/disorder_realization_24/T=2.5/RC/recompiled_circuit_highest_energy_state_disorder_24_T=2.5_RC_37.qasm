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
rz(-0.1877187) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(4.6133572) q[1];
sx q[1];
rz(7.5724966) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780802) q[0];
sx q[0];
rz(-2.8239282) q[0];
sx q[0];
rz(2.4207522) q[0];
rz(-0.046229049) q[2];
sx q[2];
rz(-2.5704489) q[2];
sx q[2];
rz(2.2276304) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6481406) q[1];
sx q[1];
rz(-1.3228184) q[1];
sx q[1];
rz(-0.47391717) q[1];
rz(-pi) q[2];
rz(-0.23203316) q[3];
sx q[3];
rz(-1.9153144) q[3];
sx q[3];
rz(-1.0330878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(2.3647986) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60748196) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(0.99824655) q[0];
rz(2.7718995) q[1];
sx q[1];
rz(-2.1001215) q[1];
sx q[1];
rz(-1.1179771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3415754) q[0];
sx q[0];
rz(-1.9722) q[0];
sx q[0];
rz(2.163351) q[0];
rz(-pi) q[1];
rz(-0.67310272) q[2];
sx q[2];
rz(-0.80922283) q[2];
sx q[2];
rz(-0.26462091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0740023) q[1];
sx q[1];
rz(-0.47110456) q[1];
sx q[1];
rz(1.8995238) q[1];
rz(-pi) q[2];
rz(-2.2023674) q[3];
sx q[3];
rz(-1.7895921) q[3];
sx q[3];
rz(-2.95524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86576858) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(-2.7562874) q[2];
rz(3.1028808) q[3];
sx q[3];
rz(-0.35274115) q[3];
sx q[3];
rz(-2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63967079) q[0];
sx q[0];
rz(-0.47698912) q[0];
sx q[0];
rz(-1.3013526) q[0];
rz(3.0838857) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(1.8147963) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0327137) q[0];
sx q[0];
rz(-2.0427867) q[0];
sx q[0];
rz(0.95979877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71452629) q[2];
sx q[2];
rz(-2.9305612) q[2];
sx q[2];
rz(-2.2843366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51242764) q[1];
sx q[1];
rz(-1.8574134) q[1];
sx q[1];
rz(0.54835876) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9671828) q[3];
sx q[3];
rz(-0.76226888) q[3];
sx q[3];
rz(1.9844733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79746276) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(2.3878035) q[2];
rz(-1.3695184) q[3];
sx q[3];
rz(-1.4621719) q[3];
sx q[3];
rz(2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(-1.3947067) q[0];
rz(-1.9649547) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-0.36697695) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0622371) q[0];
sx q[0];
rz(-1.6274165) q[0];
sx q[0];
rz(1.5039773) q[0];
x q[1];
rz(2.4796333) q[2];
sx q[2];
rz(-0.41482224) q[2];
sx q[2];
rz(2.7560459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48880348) q[1];
sx q[1];
rz(-1.3769883) q[1];
sx q[1];
rz(-2.778462) q[1];
x q[2];
rz(0.36716299) q[3];
sx q[3];
rz(-1.6510186) q[3];
sx q[3];
rz(2.1047161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7299812) q[2];
sx q[2];
rz(-1.2282164) q[2];
sx q[2];
rz(-1.8709987) q[2];
rz(-0.96108428) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(0.050475033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56080317) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(-1.8286937) q[0];
rz(1.987223) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-0.55783522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8141906) q[0];
sx q[0];
rz(-2.0923776) q[0];
sx q[0];
rz(2.5445698) q[0];
rz(-pi) q[1];
x q[1];
rz(7/(5*pi)) q[2];
sx q[2];
rz(-0.17765309) q[2];
sx q[2];
rz(-0.38933094) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.085943) q[1];
sx q[1];
rz(-2.331877) q[1];
sx q[1];
rz(1.147406) q[1];
rz(-0.46164234) q[3];
sx q[3];
rz(-1.7064693) q[3];
sx q[3];
rz(1.5380579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.64342) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(-0.15527209) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7840541) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(-2.6724755) q[0];
rz(-2.6672089) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(-2.3862086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82925382) q[0];
sx q[0];
rz(-0.27129284) q[0];
sx q[0];
rz(-1.5543227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7003661) q[2];
sx q[2];
rz(-0.6775113) q[2];
sx q[2];
rz(-0.73523075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0363706) q[1];
sx q[1];
rz(-1.3579988) q[1];
sx q[1];
rz(-0.12963055) q[1];
rz(-pi) q[2];
rz(-2.3085318) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(-1.2556094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(2.9537436) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.93194) q[0];
sx q[0];
rz(-1.2670452) q[0];
sx q[0];
rz(-0.39091045) q[0];
rz(-2.2115425) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(-1.3040868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1159853) q[0];
sx q[0];
rz(-1.5497594) q[0];
sx q[0];
rz(-1.2198835) q[0];
rz(1.6565336) q[2];
sx q[2];
rz(-1.2303599) q[2];
sx q[2];
rz(-0.96735937) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.89829373) q[1];
sx q[1];
rz(-2.0212681) q[1];
sx q[1];
rz(-0.18646481) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.061873925) q[3];
sx q[3];
rz(-0.29682595) q[3];
sx q[3];
rz(0.59943953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2685252) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(-3.0282057) q[2];
rz(-3.125627) q[3];
sx q[3];
rz(-1.4957875) q[3];
sx q[3];
rz(0.16460831) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575344) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(3.0175324) q[0];
rz(-0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(-1.442499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589813) q[0];
sx q[0];
rz(-1.9537203) q[0];
sx q[0];
rz(3.0331066) q[0];
x q[1];
rz(0.29192544) q[2];
sx q[2];
rz(-1.6124083) q[2];
sx q[2];
rz(2.6440563) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1349843) q[1];
sx q[1];
rz(-1.4378752) q[1];
sx q[1];
rz(-1.9167109) q[1];
rz(-pi) q[2];
rz(-2.8644531) q[3];
sx q[3];
rz(-2.3921514) q[3];
sx q[3];
rz(0.14271862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1759935) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(-0.74481258) q[2];
rz(0.92721573) q[3];
sx q[3];
rz(-1.6330481) q[3];
sx q[3];
rz(1.0923227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-0.48831478) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-0.22536817) q[0];
rz(2.0291746) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(0.89881277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5170636) q[0];
sx q[0];
rz(-1.0072021) q[0];
sx q[0];
rz(2.5007912) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5139047) q[2];
sx q[2];
rz(-1.2067087) q[2];
sx q[2];
rz(-0.68825005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4087993) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(0.2562457) q[1];
rz(-pi) q[2];
rz(-3.131116) q[3];
sx q[3];
rz(-1.2779405) q[3];
sx q[3];
rz(1.7550857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.044346873) q[2];
sx q[2];
rz(-1.711742) q[2];
sx q[2];
rz(-2.782235) q[2];
rz(0.25092956) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(-0.76771626) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3971685) q[0];
sx q[0];
rz(-0.1305307) q[0];
sx q[0];
rz(2.2644444) q[0];
rz(1.5646704) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(0.78561479) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79646677) q[0];
sx q[0];
rz(-0.32497999) q[0];
sx q[0];
rz(-0.82635211) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9989817) q[2];
sx q[2];
rz(-1.0773563) q[2];
sx q[2];
rz(2.7827415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0500846) q[1];
sx q[1];
rz(-0.35165241) q[1];
sx q[1];
rz(2.9356454) q[1];
rz(-pi) q[2];
rz(-1.7212058) q[3];
sx q[3];
rz(-1.0690983) q[3];
sx q[3];
rz(0.13817638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(-0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.8332558) q[3];
sx q[3];
rz(-2.8452528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(2.234266) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(0.39225929) q[2];
sx q[2];
rz(-2.1156559) q[2];
sx q[2];
rz(1.0012913) q[2];
rz(-1.0492269) q[3];
sx q[3];
rz(-2.0502301) q[3];
sx q[3];
rz(-0.32296317) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
