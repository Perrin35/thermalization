OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.719425) q[0];
sx q[0];
rz(-0.72501215) q[0];
sx q[0];
rz(-2.3056735) q[0];
rz(0.82756502) q[1];
sx q[1];
rz(-2.8711072) q[1];
sx q[1];
rz(-0.40721133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93210685) q[0];
sx q[0];
rz(-1.7833685) q[0];
sx q[0];
rz(-1.2680156) q[0];
rz(-1.9303481) q[2];
sx q[2];
rz(-1.7952732) q[2];
sx q[2];
rz(-1.4663855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8460257) q[1];
sx q[1];
rz(-0.4163308) q[1];
sx q[1];
rz(-0.33748547) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5312649) q[3];
sx q[3];
rz(-0.84926987) q[3];
sx q[3];
rz(1.5278097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.779458) q[2];
sx q[2];
rz(-0.3388277) q[2];
sx q[2];
rz(2.1753878) q[2];
rz(0.90159121) q[3];
sx q[3];
rz(-2.2880771) q[3];
sx q[3];
rz(1.0454267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40039429) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(-0.10261593) q[0];
rz(1.1688894) q[1];
sx q[1];
rz(-0.97016197) q[1];
sx q[1];
rz(2.6279367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9285787) q[0];
sx q[0];
rz(-1.4800819) q[0];
sx q[0];
rz(2.0558393) q[0];
x q[1];
rz(-1.7134708) q[2];
sx q[2];
rz(-2.1032984) q[2];
sx q[2];
rz(2.9994158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2741435) q[1];
sx q[1];
rz(-1.2701057) q[1];
sx q[1];
rz(-0.97598981) q[1];
rz(-2.0046141) q[3];
sx q[3];
rz(-1.2946715) q[3];
sx q[3];
rz(-1.9406589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1932842) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(-2.9330758) q[2];
rz(-0.68566132) q[3];
sx q[3];
rz(-1.0582358) q[3];
sx q[3];
rz(2.9764777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.0399703) q[1];
sx q[1];
rz(-1.5922348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4195815) q[0];
sx q[0];
rz(-1.7529468) q[0];
sx q[0];
rz(2.8246882) q[0];
rz(-2.9804055) q[2];
sx q[2];
rz(-1.6897795) q[2];
sx q[2];
rz(-0.42725249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0682837) q[1];
sx q[1];
rz(-1.4335634) q[1];
sx q[1];
rz(-2.2925966) q[1];
rz(-pi) q[2];
rz(-2.3089448) q[3];
sx q[3];
rz(-0.73048985) q[3];
sx q[3];
rz(-0.60527847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0478829) q[2];
sx q[2];
rz(-1.9423395) q[2];
sx q[2];
rz(-1.017978) q[2];
rz(1.7300946) q[3];
sx q[3];
rz(-0.74677765) q[3];
sx q[3];
rz(-1.5311034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3282851) q[0];
sx q[0];
rz(-1.3628553) q[0];
sx q[0];
rz(0.12411975) q[0];
rz(-0.92503754) q[1];
sx q[1];
rz(-1.8412453) q[1];
sx q[1];
rz(-0.73840028) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772241) q[0];
sx q[0];
rz(-1.508655) q[0];
sx q[0];
rz(0.81838276) q[0];
x q[1];
rz(-1.1843119) q[2];
sx q[2];
rz(-0.99683657) q[2];
sx q[2];
rz(-1.5939764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70917801) q[1];
sx q[1];
rz(-1.3644618) q[1];
sx q[1];
rz(0.56992759) q[1];
rz(-pi) q[2];
rz(-0.57598007) q[3];
sx q[3];
rz(-1.4864731) q[3];
sx q[3];
rz(0.67614851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.063252123) q[2];
sx q[2];
rz(-1.9535148) q[2];
sx q[2];
rz(-1.1032907) q[2];
rz(0.25857806) q[3];
sx q[3];
rz(-2.0662112) q[3];
sx q[3];
rz(-0.28178373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89541188) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(-1.8038764) q[0];
rz(0.61996639) q[1];
sx q[1];
rz(-1.4620616) q[1];
sx q[1];
rz(1.6372797) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02649826) q[0];
sx q[0];
rz(-2.8015602) q[0];
sx q[0];
rz(-2.1676201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4392534) q[2];
sx q[2];
rz(-1.5389256) q[2];
sx q[2];
rz(1.396871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41188365) q[1];
sx q[1];
rz(-2.4296654) q[1];
sx q[1];
rz(-0.61213778) q[1];
rz(-2.3035731) q[3];
sx q[3];
rz(-1.4921643) q[3];
sx q[3];
rz(1.3669612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6164246) q[2];
sx q[2];
rz(-1.1526356) q[2];
sx q[2];
rz(-2.0174513) q[2];
rz(0.21640402) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(-2.0719349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44806099) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(3.0060153) q[1];
sx q[1];
rz(-2.1453073) q[1];
sx q[1];
rz(2.7574976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.102589) q[0];
sx q[0];
rz(-1.6854992) q[0];
sx q[0];
rz(-2.6713598) q[0];
rz(-pi) q[1];
rz(2.9303418) q[2];
sx q[2];
rz(-2.3659035) q[2];
sx q[2];
rz(2.3084909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1138224) q[1];
sx q[1];
rz(-2.7384187) q[1];
sx q[1];
rz(1.1330963) q[1];
rz(-2.1190507) q[3];
sx q[3];
rz(-1.5060079) q[3];
sx q[3];
rz(-0.2906876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7163081) q[2];
sx q[2];
rz(-2.7649438) q[2];
sx q[2];
rz(-0.46871218) q[2];
rz(-2.5675755) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(-0.65042692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531517) q[0];
sx q[0];
rz(-1.783239) q[0];
sx q[0];
rz(0.69001946) q[0];
rz(0.72235876) q[1];
sx q[1];
rz(-2.0004309) q[1];
sx q[1];
rz(2.3648327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070097797) q[0];
sx q[0];
rz(-2.8970708) q[0];
sx q[0];
rz(-0.31425176) q[0];
rz(-pi) q[1];
rz(2.6138814) q[2];
sx q[2];
rz(-2.4832758) q[2];
sx q[2];
rz(1.1084313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0852564) q[1];
sx q[1];
rz(-1.3099226) q[1];
sx q[1];
rz(-0.44020997) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89069907) q[3];
sx q[3];
rz(-1.5560702) q[3];
sx q[3];
rz(0.17531987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7338099) q[2];
sx q[2];
rz(-1.5747384) q[2];
sx q[2];
rz(0.053827914) q[2];
rz(-2.5194061) q[3];
sx q[3];
rz(-2.6572808) q[3];
sx q[3];
rz(0.079306451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4902041) q[0];
sx q[0];
rz(-1.1505928) q[0];
sx q[0];
rz(1.401061) q[0];
rz(0.98848629) q[1];
sx q[1];
rz(-1.2683615) q[1];
sx q[1];
rz(-0.40774694) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893129) q[0];
sx q[0];
rz(-1.3527217) q[0];
sx q[0];
rz(-0.41357354) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2456535) q[2];
sx q[2];
rz(-1.6701588) q[2];
sx q[2];
rz(-2.4473913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11800471) q[1];
sx q[1];
rz(-2.0366001) q[1];
sx q[1];
rz(-0.76264571) q[1];
rz(-pi) q[2];
rz(-1.4325822) q[3];
sx q[3];
rz(-1.8402037) q[3];
sx q[3];
rz(2.9215653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72655788) q[2];
sx q[2];
rz(-2.0099535) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938195) q[0];
sx q[0];
rz(-1.6827826) q[0];
sx q[0];
rz(2.7517125) q[0];
rz(1.9436721) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(0.87497154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724204) q[0];
sx q[0];
rz(-1.5932788) q[0];
sx q[0];
rz(-2.1940986) q[0];
rz(-pi) q[1];
rz(2.2570043) q[2];
sx q[2];
rz(-1.8914127) q[2];
sx q[2];
rz(-2.5845762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5510721) q[1];
sx q[1];
rz(-1.1391907) q[1];
sx q[1];
rz(-1.1882395) q[1];
x q[2];
rz(-0.95621527) q[3];
sx q[3];
rz(-0.66171911) q[3];
sx q[3];
rz(-0.9635409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3206869) q[2];
sx q[2];
rz(-2.4497538) q[2];
sx q[2];
rz(-0.64986491) q[2];
rz(-1.1167022) q[3];
sx q[3];
rz(-1.8442804) q[3];
sx q[3];
rz(-2.9445649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7277302) q[0];
sx q[0];
rz(-0.05302269) q[0];
sx q[0];
rz(-2.9470288) q[0];
rz(-2.3692756) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(0.17759855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5935737) q[0];
sx q[0];
rz(-1.5724036) q[0];
sx q[0];
rz(-3.1409825) q[0];
rz(-pi) q[1];
rz(-0.76552154) q[2];
sx q[2];
rz(-2.3280848) q[2];
sx q[2];
rz(-3.0610656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1298123) q[1];
sx q[1];
rz(-2.3200071) q[1];
sx q[1];
rz(-0.66844861) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4836618) q[3];
sx q[3];
rz(-2.2773491) q[3];
sx q[3];
rz(2.547675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3557768) q[2];
sx q[2];
rz(-0.59712258) q[2];
sx q[2];
rz(-0.43006483) q[2];
rz(-2.0307342) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(-0.41657579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149902) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(1.968374) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(-1.0991312) q[2];
sx q[2];
rz(-0.57715125) q[2];
sx q[2];
rz(-1.1751529) q[2];
rz(-2.1803754) q[3];
sx q[3];
rz(-0.81351316) q[3];
sx q[3];
rz(-1.0143435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
