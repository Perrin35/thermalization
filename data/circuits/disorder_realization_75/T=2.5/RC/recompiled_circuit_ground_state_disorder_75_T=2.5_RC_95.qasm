OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2158382) q[0];
sx q[0];
rz(-2.821142) q[0];
sx q[0];
rz(-3.0132063) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(-0.58712062) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059537236) q[0];
sx q[0];
rz(-2.035241) q[0];
sx q[0];
rz(-1.4091253) q[0];
rz(-pi) q[1];
rz(0.70694114) q[2];
sx q[2];
rz(-0.64919186) q[2];
sx q[2];
rz(0.45237088) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93361357) q[1];
sx q[1];
rz(-0.99665239) q[1];
sx q[1];
rz(-1.0453392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5641603) q[3];
sx q[3];
rz(-0.38449461) q[3];
sx q[3];
rz(0.16669434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2653653) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(1.5246576) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(-0.18572148) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(-1.349378) q[0];
rz(-2.7941864) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.4030392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058075) q[0];
sx q[0];
rz(-0.62816915) q[0];
sx q[0];
rz(-2.7676393) q[0];
rz(-pi) q[1];
rz(-2.5591056) q[2];
sx q[2];
rz(-0.91160027) q[2];
sx q[2];
rz(2.2381353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0885873) q[1];
sx q[1];
rz(-2.2777875) q[1];
sx q[1];
rz(-1.7729575) q[1];
rz(-pi) q[2];
rz(-0.88256695) q[3];
sx q[3];
rz(-1.7808796) q[3];
sx q[3];
rz(2.9391367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8189524) q[2];
sx q[2];
rz(-1.8337269) q[2];
sx q[2];
rz(2.9827706) q[2];
rz(-2.9239376) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(-1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36667103) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(-1.2694673) q[0];
rz(-0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8458662) q[0];
sx q[0];
rz(-0.55385607) q[0];
sx q[0];
rz(-0.94828301) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4478618) q[2];
sx q[2];
rz(-0.28273928) q[2];
sx q[2];
rz(3.0764584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.060614) q[1];
sx q[1];
rz(-1.9658057) q[1];
sx q[1];
rz(-3.0837713) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9004542) q[3];
sx q[3];
rz(-1.0287549) q[3];
sx q[3];
rz(-2.295376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.082108214) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(-0.19169894) q[2];
rz(-0.55255237) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(2.1171169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315652) q[0];
sx q[0];
rz(-0.64842328) q[0];
sx q[0];
rz(-0.60436526) q[0];
rz(-1.2454237) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(0.85974685) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5167546) q[0];
sx q[0];
rz(-2.4167912) q[0];
sx q[0];
rz(-1.047931) q[0];
rz(-1.672869) q[2];
sx q[2];
rz(-1.5849176) q[2];
sx q[2];
rz(0.532224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8670876) q[1];
sx q[1];
rz(-0.18840677) q[1];
sx q[1];
rz(2.8727358) q[1];
rz(0.28713496) q[3];
sx q[3];
rz(-2.4004712) q[3];
sx q[3];
rz(2.5484619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2802281) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0786521) q[0];
sx q[0];
rz(-0.27314726) q[0];
sx q[0];
rz(2.8237421) q[0];
rz(-1.3392797) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(0.97871614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1237549) q[0];
sx q[0];
rz(-0.40214254) q[0];
sx q[0];
rz(3.0642302) q[0];
rz(1.481858) q[2];
sx q[2];
rz(-1.3761889) q[2];
sx q[2];
rz(-0.22502514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3826344) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(-1.7428223) q[1];
rz(-pi) q[2];
rz(2.948165) q[3];
sx q[3];
rz(-1.7599835) q[3];
sx q[3];
rz(2.6003169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19834441) q[2];
sx q[2];
rz(-0.94503108) q[2];
sx q[2];
rz(1.1172969) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(1.1972637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4365874) q[0];
sx q[0];
rz(-2.8745108) q[0];
sx q[0];
rz(1.3096814) q[0];
rz(1.0218703) q[1];
sx q[1];
rz(-0.836687) q[1];
sx q[1];
rz(-1.4885363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90327901) q[0];
sx q[0];
rz(-0.92172232) q[0];
sx q[0];
rz(0.42239503) q[0];
x q[1];
rz(2.4632881) q[2];
sx q[2];
rz(-1.3875752) q[2];
sx q[2];
rz(-1.2228633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0206581) q[1];
sx q[1];
rz(-1.7132732) q[1];
sx q[1];
rz(-2.6441036) q[1];
rz(-pi) q[2];
x q[2];
rz(1.009935) q[3];
sx q[3];
rz(-1.4743685) q[3];
sx q[3];
rz(1.7280462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66270193) q[2];
sx q[2];
rz(-1.3871437) q[2];
sx q[2];
rz(0.26408163) q[2];
rz(-1.4763907) q[3];
sx q[3];
rz(-0.68015209) q[3];
sx q[3];
rz(-2.7217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61178094) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(1.06426) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(1.4605716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34188627) q[0];
sx q[0];
rz(-2.0545022) q[0];
sx q[0];
rz(1.6428436) q[0];
rz(-pi) q[1];
rz(1.2719897) q[2];
sx q[2];
rz(-1.5332216) q[2];
sx q[2];
rz(0.36848289) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8014508) q[1];
sx q[1];
rz(-0.51348084) q[1];
sx q[1];
rz(2.5380847) q[1];
rz(2.0123294) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(-0.03015524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(2.5420945) q[2];
rz(0.89764578) q[3];
sx q[3];
rz(-1.7220327) q[3];
sx q[3];
rz(-3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(-1.7401485) q[0];
rz(-3.0701045) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(2.1404526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6342696) q[0];
sx q[0];
rz(-1.3183013) q[0];
sx q[0];
rz(0.4613614) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1835353) q[2];
sx q[2];
rz(-0.94410482) q[2];
sx q[2];
rz(-1.3037217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65232507) q[1];
sx q[1];
rz(-1.7026002) q[1];
sx q[1];
rz(2.5062923) q[1];
rz(-0.66949943) q[3];
sx q[3];
rz(-1.0855128) q[3];
sx q[3];
rz(-0.2822596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1772168) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(3.0948203) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(-2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319594) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(0.69044789) q[0];
rz(-0.18028232) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(2.8188425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0092516) q[0];
sx q[0];
rz(-1.7242431) q[0];
sx q[0];
rz(-1.4645534) q[0];
rz(-pi) q[1];
rz(1.5602925) q[2];
sx q[2];
rz(-0.5521419) q[2];
sx q[2];
rz(-0.10607468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.018498713) q[1];
sx q[1];
rz(-1.0375751) q[1];
sx q[1];
rz(-0.64792222) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2261103) q[3];
sx q[3];
rz(-0.99952664) q[3];
sx q[3];
rz(-0.20184982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6705769) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(-2.1237109) q[2];
rz(0.0035303591) q[3];
sx q[3];
rz(-2.1719833) q[3];
sx q[3];
rz(1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8135391) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(0.92593431) q[0];
rz(-3.0042341) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(2.5416809) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900676) q[0];
sx q[0];
rz(-2.5976564) q[0];
sx q[0];
rz(-0.005225709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0157064) q[2];
sx q[2];
rz(-1.0169463) q[2];
sx q[2];
rz(-0.24283838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.74739) q[1];
sx q[1];
rz(-2.5521186) q[1];
sx q[1];
rz(1.6255476) q[1];
rz(-pi) q[2];
rz(0.32528721) q[3];
sx q[3];
rz(-1.1184177) q[3];
sx q[3];
rz(-0.77879209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4124734) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(-0.64175433) q[2];
rz(1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(-1.6183287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625576) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(1.7017801) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-2.4859602) q[2];
sx q[2];
rz(-0.31972319) q[2];
sx q[2];
rz(2.2421851) q[2];
rz(2.6445848) q[3];
sx q[3];
rz(-1.0773226) q[3];
sx q[3];
rz(1.1763431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
