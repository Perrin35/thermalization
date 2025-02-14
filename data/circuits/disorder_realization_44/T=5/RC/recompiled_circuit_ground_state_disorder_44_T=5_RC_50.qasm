OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(-2.1679572) q[1];
sx q[1];
rz(-1.4373625) q[1];
sx q[1];
rz(0.91926423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96199233) q[0];
sx q[0];
rz(-1.5536947) q[0];
sx q[0];
rz(0.03058612) q[0];
rz(-pi) q[1];
rz(-0.12030258) q[2];
sx q[2];
rz(-1.6610816) q[2];
sx q[2];
rz(1.1995969) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8074933) q[1];
sx q[1];
rz(-1.6624644) q[1];
sx q[1];
rz(-0.86218545) q[1];
rz(-pi) q[2];
rz(2.0415061) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(1.4130886) q[2];
rz(-0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(0.10281674) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8973812) q[0];
sx q[0];
rz(-2.057071) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(2.3731025) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(-2.1220727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1508038) q[0];
sx q[0];
rz(-1.6103364) q[0];
sx q[0];
rz(-2.6943452) q[0];
x q[1];
rz(-2.0020194) q[2];
sx q[2];
rz(-1.5016455) q[2];
sx q[2];
rz(0.56996843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63377127) q[1];
sx q[1];
rz(-1.941676) q[1];
sx q[1];
rz(0.033286496) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7048265) q[3];
sx q[3];
rz(-0.55826) q[3];
sx q[3];
rz(-0.2414862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(-1.2223318) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(1.2179751) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(0.9224433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24127125) q[0];
sx q[0];
rz(-1.5301955) q[0];
sx q[0];
rz(-3.0837644) q[0];
x q[1];
rz(-1.0664682) q[2];
sx q[2];
rz(-2.2309003) q[2];
sx q[2];
rz(2.4911027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53239142) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(0.81873853) q[1];
rz(-pi) q[2];
rz(0.32914583) q[3];
sx q[3];
rz(-1.4107879) q[3];
sx q[3];
rz(0.015586675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.018365232) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(1.3398735) q[2];
rz(2.5943622) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0860586) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-0.25746447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5437357) q[0];
sx q[0];
rz(-0.86559767) q[0];
sx q[0];
rz(-1.0969093) q[0];
rz(-pi) q[1];
rz(-0.065071062) q[2];
sx q[2];
rz(-2.3461968) q[2];
sx q[2];
rz(-0.8000904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83492571) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(2.5132781) q[1];
x q[2];
rz(-2.5268528) q[3];
sx q[3];
rz(-2.2343544) q[3];
sx q[3];
rz(-3.1082982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(0.47284687) q[2];
rz(-2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351629) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(3.0761062) q[0];
rz(0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.4261036) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55520362) q[0];
sx q[0];
rz(-2.9034333) q[0];
sx q[0];
rz(2.737153) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21435301) q[2];
sx q[2];
rz(-2.7446236) q[2];
sx q[2];
rz(3.1052542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8219442) q[1];
sx q[1];
rz(-1.8488171) q[1];
sx q[1];
rz(1.3389498) q[1];
rz(0.17284837) q[3];
sx q[3];
rz(-1.2842872) q[3];
sx q[3];
rz(0.89927538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19491974) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(1.8348414) q[2];
rz(-1.170018) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829247) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(-0.40147993) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-0.79877001) q[1];
sx q[1];
rz(-2.2875517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14388785) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(-1.3530988) q[0];
rz(-2.4638205) q[2];
sx q[2];
rz(-0.23636625) q[2];
sx q[2];
rz(1.4217699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6352977) q[1];
sx q[1];
rz(-0.98165252) q[1];
sx q[1];
rz(-1.4511257) q[1];
rz(1.1229418) q[3];
sx q[3];
rz(-1.0077701) q[3];
sx q[3];
rz(-0.26859586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4868769) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(1.9640478) q[2];
rz(-2.5901637) q[3];
sx q[3];
rz(-1.5588372) q[3];
sx q[3];
rz(2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42171445) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-0.91845671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699444) q[0];
sx q[0];
rz(-1.3027096) q[0];
sx q[0];
rz(1.2275342) q[0];
x q[1];
rz(2.2662042) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(0.018176807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.368093) q[1];
sx q[1];
rz(-1.5640904) q[1];
sx q[1];
rz(-3.1342602) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9711668) q[3];
sx q[3];
rz(-1.2316576) q[3];
sx q[3];
rz(-2.2234774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.034417001) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(1.3431965) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69304943) q[0];
sx q[0];
rz(-2.9149084) q[0];
sx q[0];
rz(1.6402798) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728031) q[0];
sx q[0];
rz(-0.69376341) q[0];
sx q[0];
rz(-2.540178) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2605687) q[2];
sx q[2];
rz(-2.0379279) q[2];
sx q[2];
rz(1.9595944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4290165) q[1];
sx q[1];
rz(-1.4498931) q[1];
sx q[1];
rz(2.4790133) q[1];
rz(-0.39466484) q[3];
sx q[3];
rz(-1.695444) q[3];
sx q[3];
rz(0.11388557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1559653) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(-2.8333832) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(-1.6109899) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2823328) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(-1.169091) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(-1.4287359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7633782) q[0];
sx q[0];
rz(-1.175048) q[0];
sx q[0];
rz(-1.2984492) q[0];
rz(2.4506017) q[2];
sx q[2];
rz(-0.96110247) q[2];
sx q[2];
rz(0.26870773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1395806) q[1];
sx q[1];
rz(-1.5095469) q[1];
sx q[1];
rz(-1.9467926) q[1];
rz(0.47023021) q[3];
sx q[3];
rz(-2.2080126) q[3];
sx q[3];
rz(-1.5772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94065654) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(-1.3915871) q[2];
rz(0.40677795) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(-0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(-1.9039924) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(-0.79992574) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7800723) q[0];
sx q[0];
rz(-1.8006386) q[0];
sx q[0];
rz(2.337268) q[0];
rz(-pi) q[1];
rz(2.6705856) q[2];
sx q[2];
rz(-2.5386497) q[2];
sx q[2];
rz(-2.431536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7042735) q[1];
sx q[1];
rz(-1.4156439) q[1];
sx q[1];
rz(-0.67351933) q[1];
rz(-2.0560451) q[3];
sx q[3];
rz(-1.3486514) q[3];
sx q[3];
rz(-2.9911161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7633729) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(0.073089449) q[2];
rz(-2.8857005) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(-1.5727829) q[2];
sx q[2];
rz(-2.28021) q[2];
sx q[2];
rz(1.8328666) q[2];
rz(-0.70684915) q[3];
sx q[3];
rz(-0.77597386) q[3];
sx q[3];
rz(1.5472277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
