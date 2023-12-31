OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(-2.4175329) q[0];
sx q[0];
rz(1.6568503) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523525) q[0];
sx q[0];
rz(-1.7416735) q[0];
sx q[0];
rz(1.6160374) q[0];
rz(-pi) q[1];
rz(-2.9759679) q[2];
sx q[2];
rz(-2.1016444) q[2];
sx q[2];
rz(2.3373375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29477316) q[1];
sx q[1];
rz(-1.5713437) q[1];
sx q[1];
rz(0.16695395) q[1];
rz(-pi) q[2];
rz(3.1021318) q[3];
sx q[3];
rz(-0.69898116) q[3];
sx q[3];
rz(0.82074245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(-3.0956691) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(-1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-3.0156946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4100285) q[0];
sx q[0];
rz(-1.0985939) q[0];
sx q[0];
rz(-2.7239951) q[0];
rz(-2.5198031) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(0.34678005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0071348) q[1];
sx q[1];
rz(-1.6975228) q[1];
sx q[1];
rz(0.030220672) q[1];
x q[2];
rz(1.5829854) q[3];
sx q[3];
rz(-0.72353957) q[3];
sx q[3];
rz(-2.7140868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(2.753567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(-0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(-3.0575867) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(-1.9817339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1201598) q[0];
sx q[0];
rz(-0.5172356) q[0];
sx q[0];
rz(-1.4034127) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6339176) q[2];
sx q[2];
rz(-1.0055563) q[2];
sx q[2];
rz(-0.6914247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0594271) q[1];
sx q[1];
rz(-1.1929999) q[1];
sx q[1];
rz(2.6487745) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0852774) q[3];
sx q[3];
rz(-2.115613) q[3];
sx q[3];
rz(-0.69308263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(0.28465095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68344224) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(0.083806888) q[0];
rz(-0.95732032) q[2];
sx q[2];
rz(-0.34733221) q[2];
sx q[2];
rz(0.50328244) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.804867) q[1];
sx q[1];
rz(-1.4443828) q[1];
sx q[1];
rz(0.079041914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90162006) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(-1.3265142) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(-0.93793905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4612761) q[0];
sx q[0];
rz(-3.0468867) q[0];
sx q[0];
rz(1.8041496) q[0];
rz(0.51860923) q[2];
sx q[2];
rz(-1.0728288) q[2];
sx q[2];
rz(-1.2155611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.282498) q[1];
sx q[1];
rz(-1.5364093) q[1];
sx q[1];
rz(-3.1153468) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35477562) q[3];
sx q[3];
rz(-1.4491023) q[3];
sx q[3];
rz(-1.0697168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43907169) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(-2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(2.5851137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300886) q[0];
sx q[0];
rz(-0.43404365) q[0];
sx q[0];
rz(-2.6326798) q[0];
rz(0.85530497) q[2];
sx q[2];
rz(-1.7024634) q[2];
sx q[2];
rz(1.635074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1408128) q[1];
sx q[1];
rz(-1.3197348) q[1];
sx q[1];
rz(2.954133) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3658386) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(-3.0690103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(1.2247941) q[2];
rz(1.6312284) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(2.6409805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4780873) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(1.7630793) q[0];
rz(-2.0666276) q[2];
sx q[2];
rz(-0.60056409) q[2];
sx q[2];
rz(1.9753089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(1.7016181) q[1];
rz(1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(1.0637103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.0764517) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(1.4656461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38795234) q[0];
sx q[0];
rz(-2.6590829) q[0];
sx q[0];
rz(1.2778736) q[0];
x q[1];
rz(-1.4618116) q[2];
sx q[2];
rz(-1.3853067) q[2];
sx q[2];
rz(0.89599228) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17897478) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(2.7589873) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-1.0675758) q[3];
sx q[3];
rz(-2.7748231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(-2.6749558) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-0.73721686) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(2.2163056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7172456) q[0];
sx q[0];
rz(-2.2204917) q[0];
sx q[0];
rz(-0.93066494) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66321744) q[2];
sx q[2];
rz(-1.6753917) q[2];
sx q[2];
rz(1.6246206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8268938) q[1];
sx q[1];
rz(-1.3876545) q[1];
sx q[1];
rz(-1.7032743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0989283) q[3];
sx q[3];
rz(-2.3657551) q[3];
sx q[3];
rz(0.27234205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(-1.0037237) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(0.75751799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.950338) q[0];
sx q[0];
rz(-1.3779103) q[0];
sx q[0];
rz(-0.14001503) q[0];
rz(2.9847758) q[2];
sx q[2];
rz(-0.98625253) q[2];
sx q[2];
rz(-2.2371694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65749189) q[1];
sx q[1];
rz(-1.3490281) q[1];
sx q[1];
rz(1.1360672) q[1];
rz(3.1086139) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(-2.9041293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929852) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(-2.7040504) q[2];
sx q[2];
rz(-1.2699288) q[2];
sx q[2];
rz(-1.5581836) q[2];
rz(-2.8201841) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
