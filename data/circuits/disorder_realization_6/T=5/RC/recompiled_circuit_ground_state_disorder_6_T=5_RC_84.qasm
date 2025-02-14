OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(0.21831231) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(-0.47774878) q[1];
sx q[1];
rz(-2.5876317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4628578) q[0];
sx q[0];
rz(-2.4241872) q[0];
sx q[0];
rz(-0.82490246) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6126819) q[2];
sx q[2];
rz(-1.5195091) q[2];
sx q[2];
rz(-1.9330658) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40383717) q[1];
sx q[1];
rz(-0.74450508) q[1];
sx q[1];
rz(-1.2584524) q[1];
rz(-pi) q[2];
rz(-1.9945108) q[3];
sx q[3];
rz(-1.0726352) q[3];
sx q[3];
rz(-0.26396449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7654984) q[2];
sx q[2];
rz(-2.8480397) q[2];
sx q[2];
rz(1.2676839) q[2];
rz(0.82967657) q[3];
sx q[3];
rz(-1.6206348) q[3];
sx q[3];
rz(2.7177496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.394519) q[0];
rz(2.246619) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5590394) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48838636) q[0];
sx q[0];
rz(-2.1108642) q[0];
sx q[0];
rz(-0.77629838) q[0];
x q[1];
rz(1.1505914) q[2];
sx q[2];
rz(-1.8981427) q[2];
sx q[2];
rz(-0.98132747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29702052) q[1];
sx q[1];
rz(-2.3851042) q[1];
sx q[1];
rz(-0.46539657) q[1];
x q[2];
rz(-0.96554324) q[3];
sx q[3];
rz(-2.6355834) q[3];
sx q[3];
rz(0.11270302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47722694) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(-0.030755432) q[2];
rz(0.45298806) q[3];
sx q[3];
rz(-0.22946295) q[3];
sx q[3];
rz(-0.17791137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40193108) q[0];
sx q[0];
rz(-0.10098305) q[0];
sx q[0];
rz(-2.3133551) q[0];
rz(-3.0939057) q[1];
sx q[1];
rz(-0.86367718) q[1];
sx q[1];
rz(1.9140859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9253474) q[0];
sx q[0];
rz(-2.0106533) q[0];
sx q[0];
rz(0.60103215) q[0];
x q[1];
rz(-2.0426072) q[2];
sx q[2];
rz(-2.6919439) q[2];
sx q[2];
rz(1.4020021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32900287) q[1];
sx q[1];
rz(-2.0101476) q[1];
sx q[1];
rz(-2.8529608) q[1];
x q[2];
rz(-1.5121721) q[3];
sx q[3];
rz(-0.72949648) q[3];
sx q[3];
rz(0.54704715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.092827169) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(-2.0406593) q[2];
rz(0.01384211) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(-2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58532995) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(-2.8072667) q[0];
rz(-2.4090134) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(-0.11071959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78149139) q[0];
sx q[0];
rz(-1.5859787) q[0];
sx q[0];
rz(2.380037) q[0];
rz(-2.9760095) q[2];
sx q[2];
rz(-1.6628169) q[2];
sx q[2];
rz(-1.657287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19921215) q[1];
sx q[1];
rz(-0.95006493) q[1];
sx q[1];
rz(-2.278028) q[1];
x q[2];
rz(-2.0483253) q[3];
sx q[3];
rz(-0.85460409) q[3];
sx q[3];
rz(-1.0323857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4577786) q[2];
sx q[2];
rz(-2.8595698) q[2];
sx q[2];
rz(-1.7337743) q[2];
rz(1.9862566) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69283501) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(2.1441929) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(1.4195199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888694) q[0];
sx q[0];
rz(-0.27452454) q[0];
sx q[0];
rz(-2.2122266) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77395524) q[2];
sx q[2];
rz(-1.4593235) q[2];
sx q[2];
rz(2.39448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3910663) q[1];
sx q[1];
rz(-1.2519893) q[1];
sx q[1];
rz(-0.90853779) q[1];
rz(1.6035242) q[3];
sx q[3];
rz(-1.6943036) q[3];
sx q[3];
rz(-2.0211969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2898966) q[2];
sx q[2];
rz(-3.0471314) q[2];
sx q[2];
rz(0.32290253) q[2];
rz(-1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36393976) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(2.3349578) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-2.0239425) q[1];
sx q[1];
rz(-1.9516099) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4247852) q[0];
sx q[0];
rz(-1.2233226) q[0];
sx q[0];
rz(3.0106972) q[0];
x q[1];
rz(0.34752589) q[2];
sx q[2];
rz(-0.56171562) q[2];
sx q[2];
rz(0.39081854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8998225) q[1];
sx q[1];
rz(-0.17647753) q[1];
sx q[1];
rz(1.2736257) q[1];
rz(-2.4111536) q[3];
sx q[3];
rz(-1.4439266) q[3];
sx q[3];
rz(-0.74339429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7334062) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(-2.9564986) q[2];
rz(-1.6144276) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1899034) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(-0.21251799) q[0];
rz(0.7849794) q[1];
sx q[1];
rz(-1.4947944) q[1];
sx q[1];
rz(-0.77883887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3991645) q[0];
sx q[0];
rz(-1.0013097) q[0];
sx q[0];
rz(-0.67230255) q[0];
rz(0.53063993) q[2];
sx q[2];
rz(-2.4839249) q[2];
sx q[2];
rz(-0.69537698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0743588) q[1];
sx q[1];
rz(-1.7869084) q[1];
sx q[1];
rz(-1.3098148) q[1];
rz(-1.0209342) q[3];
sx q[3];
rz(-2.3519197) q[3];
sx q[3];
rz(1.1268953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51398858) q[2];
sx q[2];
rz(-2.6476761) q[2];
sx q[2];
rz(0.40840515) q[2];
rz(2.4397395) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34847611) q[0];
sx q[0];
rz(-1.9261253) q[0];
sx q[0];
rz(0.066019639) q[0];
rz(1.6325715) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(-2.1836233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6867179) q[0];
sx q[0];
rz(-1.2871847) q[0];
sx q[0];
rz(-2.6647749) q[0];
x q[1];
rz(0.52131781) q[2];
sx q[2];
rz(-0.88721878) q[2];
sx q[2];
rz(-0.41022656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4730395) q[1];
sx q[1];
rz(-2.2277692) q[1];
sx q[1];
rz(1.1725575) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.041298) q[3];
sx q[3];
rz(-1.417932) q[3];
sx q[3];
rz(-0.66674846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3079188) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(1.4109122) q[2];
rz(2.446567) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(-3.1237349) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0787635) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(2.1516946) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(-1.90082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.038886) q[0];
sx q[0];
rz(-2.8409344) q[0];
sx q[0];
rz(-1.3749529) q[0];
rz(-2.9389589) q[2];
sx q[2];
rz(-1.1362193) q[2];
sx q[2];
rz(-3.0960577) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1009192) q[1];
sx q[1];
rz(-1.2881491) q[1];
sx q[1];
rz(1.8428749) q[1];
rz(-pi) q[2];
rz(0.95043358) q[3];
sx q[3];
rz(-1.9470805) q[3];
sx q[3];
rz(-0.84445124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3339633) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(-0.71869746) q[2];
rz(-1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-2.1271472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(0.54866791) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(-2.3396709) q[0];
rz(-1.0879263) q[1];
sx q[1];
rz(-1.8049003) q[1];
sx q[1];
rz(-0.71802872) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1205667) q[0];
sx q[0];
rz(-0.88827288) q[0];
sx q[0];
rz(1.989407) q[0];
x q[1];
rz(-2.0640578) q[2];
sx q[2];
rz(-2.0414957) q[2];
sx q[2];
rz(2.5483709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40945986) q[1];
sx q[1];
rz(-1.902305) q[1];
sx q[1];
rz(-3.0863239) q[1];
rz(-2.4122756) q[3];
sx q[3];
rz(-2.1538247) q[3];
sx q[3];
rz(-2.5860525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3512909) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(2.6289319) q[2];
rz(-2.3971108) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(-0.7640394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22513334) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(-0.71612877) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(3.0307583) q[2];
sx q[2];
rz(-1.3468942) q[2];
sx q[2];
rz(-0.46769618) q[2];
rz(1.2848787) q[3];
sx q[3];
rz(-0.34964041) q[3];
sx q[3];
rz(1.6007363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
