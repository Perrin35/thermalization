OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(-2.4165805) q[0];
sx q[0];
rz(2.3056735) q[0];
rz(-2.3140276) q[1];
sx q[1];
rz(-0.27048549) q[1];
sx q[1];
rz(-2.7343813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4370928) q[0];
sx q[0];
rz(-1.2750393) q[0];
sx q[0];
rz(-2.9192135) q[0];
rz(-pi) q[1];
rz(1.2112445) q[2];
sx q[2];
rz(-1.3463194) q[2];
sx q[2];
rz(1.4663855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5558124) q[1];
sx q[1];
rz(-1.4364873) q[1];
sx q[1];
rz(2.7463169) q[1];
rz(-pi) q[2];
rz(2.3916733) q[3];
sx q[3];
rz(-1.1259606) q[3];
sx q[3];
rz(0.47575853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.779458) q[2];
sx q[2];
rz(-2.802765) q[2];
sx q[2];
rz(2.1753878) q[2];
rz(-0.90159121) q[3];
sx q[3];
rz(-2.2880771) q[3];
sx q[3];
rz(2.0961659) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7411984) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(-0.10261593) q[0];
rz(1.9727033) q[1];
sx q[1];
rz(-2.1714307) q[1];
sx q[1];
rz(2.6279367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1875604) q[0];
sx q[0];
rz(-0.49278773) q[0];
sx q[0];
rz(1.7634747) q[0];
rz(-2.9048237) q[2];
sx q[2];
rz(-2.5920923) q[2];
sx q[2];
rz(-0.41791805) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8508965) q[1];
sx q[1];
rz(-2.4833931) q[1];
sx q[1];
rz(2.0762411) q[1];
rz(-pi) q[2];
rz(-1.1369785) q[3];
sx q[3];
rz(-1.2946715) q[3];
sx q[3];
rz(1.9406589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9483084) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(-0.2085169) q[2];
rz(-2.4559313) q[3];
sx q[3];
rz(-2.0833569) q[3];
sx q[3];
rz(-0.16511495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044428069) q[0];
sx q[0];
rz(-1.2760289) q[0];
sx q[0];
rz(-1.9097419) q[0];
rz(-0.34230289) q[1];
sx q[1];
rz(-1.1016223) q[1];
sx q[1];
rz(1.5922348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72201113) q[0];
sx q[0];
rz(-1.7529468) q[0];
sx q[0];
rz(-2.8246882) q[0];
rz(-pi) q[1];
rz(-0.16118717) q[2];
sx q[2];
rz(-1.4518132) q[2];
sx q[2];
rz(2.7143402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.07330896) q[1];
sx q[1];
rz(-1.7080293) q[1];
sx q[1];
rz(0.84899606) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5991159) q[3];
sx q[3];
rz(-1.0546135) q[3];
sx q[3];
rz(-0.27942785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0937097) q[2];
sx q[2];
rz(-1.9423395) q[2];
sx q[2];
rz(-1.017978) q[2];
rz(-1.7300946) q[3];
sx q[3];
rz(-0.74677765) q[3];
sx q[3];
rz(-1.6104893) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3282851) q[0];
sx q[0];
rz(-1.7787373) q[0];
sx q[0];
rz(-0.12411975) q[0];
rz(-0.92503754) q[1];
sx q[1];
rz(-1.8412453) q[1];
sx q[1];
rz(-0.73840028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8014073) q[0];
sx q[0];
rz(-0.75447318) q[0];
sx q[0];
rz(-1.4800001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5321711) q[2];
sx q[2];
rz(-1.2487234) q[2];
sx q[2];
rz(-2.9009478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1494245) q[1];
sx q[1];
rz(-2.1271884) q[1];
sx q[1];
rz(-1.8144649) q[1];
rz(-pi) q[2];
rz(-2.9876338) q[3];
sx q[3];
rz(-2.5601644) q[3];
sx q[3];
rz(2.1179782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0783405) q[2];
sx q[2];
rz(-1.9535148) q[2];
sx q[2];
rz(1.1032907) q[2];
rz(0.25857806) q[3];
sx q[3];
rz(-2.0662112) q[3];
sx q[3];
rz(-0.28178373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2461808) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(-1.3377162) q[0];
rz(2.5216263) q[1];
sx q[1];
rz(-1.679531) q[1];
sx q[1];
rz(1.6372797) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1150944) q[0];
sx q[0];
rz(-2.8015602) q[0];
sx q[0];
rz(-0.97397255) q[0];
rz(1.5290544) q[2];
sx q[2];
rz(-2.2727059) q[2];
sx q[2];
rz(-0.14696444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4940987) q[1];
sx q[1];
rz(-1.955619) q[1];
sx q[1];
rz(0.61487756) q[1];
rz(0.83801954) q[3];
sx q[3];
rz(-1.6494284) q[3];
sx q[3];
rz(1.7746314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6164246) q[2];
sx q[2];
rz(-1.9889571) q[2];
sx q[2];
rz(1.1241414) q[2];
rz(-0.21640402) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(-1.0696577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6935317) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(-3.0060153) q[1];
sx q[1];
rz(-0.99628535) q[1];
sx q[1];
rz(2.7574976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039003619) q[0];
sx q[0];
rz(-1.4560934) q[0];
sx q[0];
rz(-0.47023288) q[0];
rz(0.76445396) q[2];
sx q[2];
rz(-1.4234418) q[2];
sx q[2];
rz(0.58576983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0277703) q[1];
sx q[1];
rz(-0.40317391) q[1];
sx q[1];
rz(2.0084963) q[1];
rz(-0.075875326) q[3];
sx q[3];
rz(-2.11777) q[3];
sx q[3];
rz(-1.8219655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7163081) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(-2.6728805) q[2];
rz(0.57401717) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(-0.65042692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18844093) q[0];
sx q[0];
rz(-1.3583536) q[0];
sx q[0];
rz(2.4515732) q[0];
rz(0.72235876) q[1];
sx q[1];
rz(-1.1411618) q[1];
sx q[1];
rz(-2.3648327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070097797) q[0];
sx q[0];
rz(-0.24452183) q[0];
sx q[0];
rz(0.31425176) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58906196) q[2];
sx q[2];
rz(-1.8839594) q[2];
sx q[2];
rz(2.2472266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36464035) q[1];
sx q[1];
rz(-1.1464652) q[1];
sx q[1];
rz(-1.2838497) q[1];
rz(-pi) q[2];
rz(2.2508936) q[3];
sx q[3];
rz(-1.5560702) q[3];
sx q[3];
rz(-2.9662728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40778273) q[2];
sx q[2];
rz(-1.5747384) q[2];
sx q[2];
rz(3.0877647) q[2];
rz(0.6221866) q[3];
sx q[3];
rz(-2.6572808) q[3];
sx q[3];
rz(0.079306451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902041) q[0];
sx q[0];
rz(-1.9909998) q[0];
sx q[0];
rz(1.401061) q[0];
rz(2.1531064) q[1];
sx q[1];
rz(-1.2683615) q[1];
sx q[1];
rz(0.40774694) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177454) q[0];
sx q[0];
rz(-1.9739986) q[0];
sx q[0];
rz(1.8082331) q[0];
rz(-pi) q[1];
rz(-0.8959391) q[2];
sx q[2];
rz(-1.6701588) q[2];
sx q[2];
rz(2.4473913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11800471) q[1];
sx q[1];
rz(-1.1049926) q[1];
sx q[1];
rz(-0.76264571) q[1];
rz(-pi) q[2];
rz(1.4325822) q[3];
sx q[3];
rz(-1.301389) q[3];
sx q[3];
rz(2.9215653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4150348) q[2];
sx q[2];
rz(-1.1316391) q[2];
sx q[2];
rz(2.4416907) q[2];
rz(-2.0578201) q[3];
sx q[3];
rz(-0.86728573) q[3];
sx q[3];
rz(-2.9054902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938195) q[0];
sx q[0];
rz(-1.45881) q[0];
sx q[0];
rz(-2.7517125) q[0];
rz(-1.1979206) q[1];
sx q[1];
rz(-3.0471264) q[1];
sx q[1];
rz(-0.87497154) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9720249) q[0];
sx q[0];
rz(-2.1939169) q[0];
sx q[0];
rz(0.02768691) q[0];
x q[1];
rz(1.0880555) q[2];
sx q[2];
rz(-0.74626479) q[2];
sx q[2];
rz(-1.3811228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81351133) q[1];
sx q[1];
rz(-1.9167148) q[1];
sx q[1];
rz(0.4608058) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1853774) q[3];
sx q[3];
rz(-2.4798735) q[3];
sx q[3];
rz(-2.1780518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3206869) q[2];
sx q[2];
rz(-0.69183886) q[2];
sx q[2];
rz(-0.64986491) q[2];
rz(1.1167022) q[3];
sx q[3];
rz(-1.8442804) q[3];
sx q[3];
rz(-0.19702774) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4138625) q[0];
sx q[0];
rz(-0.05302269) q[0];
sx q[0];
rz(2.9470288) q[0];
rz(-0.77231705) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(2.9639941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307482) q[0];
sx q[0];
rz(-3.1398735) q[0];
sx q[0];
rz(1.2079713) q[0];
rz(2.4899876) q[2];
sx q[2];
rz(-1.043104) q[2];
sx q[2];
rz(2.0737632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87093052) q[1];
sx q[1];
rz(-2.1829603) q[1];
sx q[1];
rz(-0.98302676) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74754369) q[3];
sx q[3];
rz(-1.0870544) q[3];
sx q[3];
rz(-1.699699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3557768) q[2];
sx q[2];
rz(-0.59712258) q[2];
sx q[2];
rz(0.43006483) q[2];
rz(-2.0307342) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(-0.41657579) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4149902) q[0];
sx q[0];
rz(-1.5885329) q[0];
sx q[0];
rz(2.9581099) q[0];
rz(-1.968374) q[1];
sx q[1];
rz(-1.9216187) q[1];
sx q[1];
rz(0.12513195) q[1];
rz(2.0424615) q[2];
sx q[2];
rz(-0.57715125) q[2];
sx q[2];
rz(-1.1751529) q[2];
rz(2.1803754) q[3];
sx q[3];
rz(-2.3280795) q[3];
sx q[3];
rz(2.1272492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
