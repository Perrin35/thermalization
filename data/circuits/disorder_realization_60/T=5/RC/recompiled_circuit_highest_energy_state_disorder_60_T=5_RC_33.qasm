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
rz(-1.6835535) q[0];
sx q[0];
rz(-2.2245421) q[0];
sx q[0];
rz(-0.24669692) q[0];
rz(3.9858272) q[1];
sx q[1];
rz(1.2935473) q[1];
sx q[1];
rz(10.727439) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3299574) q[0];
sx q[0];
rz(-2.6751452) q[0];
sx q[0];
rz(0.88282013) q[0];
rz(-pi) q[1];
rz(2.1470492) q[2];
sx q[2];
rz(-2.3495397) q[2];
sx q[2];
rz(0.22135569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6410632) q[1];
sx q[1];
rz(-1.682356) q[1];
sx q[1];
rz(2.5140933) q[1];
rz(-pi) q[2];
rz(1.4250303) q[3];
sx q[3];
rz(-1.4633388) q[3];
sx q[3];
rz(-1.2821357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.777433) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(1.1945266) q[2];
rz(1.901769) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(1.4136081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164923) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-0.66705739) q[0];
rz(0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(-1.0557231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919843) q[0];
sx q[0];
rz(-2.4469564) q[0];
sx q[0];
rz(1.9978092) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.021432545) q[2];
sx q[2];
rz(-1.1864682) q[2];
sx q[2];
rz(0.0057980428) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2438587) q[1];
sx q[1];
rz(-1.4816726) q[1];
sx q[1];
rz(-1.9088163) q[1];
x q[2];
rz(2.3930644) q[3];
sx q[3];
rz(-2.1338042) q[3];
sx q[3];
rz(-0.23458086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(-2.5248027) q[2];
rz(2.8957497) q[3];
sx q[3];
rz(-1.5522233) q[3];
sx q[3];
rz(1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-0.73391947) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(-0.47897419) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(-1.3165547) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.115057) q[0];
sx q[0];
rz(-2.2000072) q[0];
sx q[0];
rz(-2.3903484) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9390305) q[2];
sx q[2];
rz(-1.1656467) q[2];
sx q[2];
rz(-0.7140401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41539792) q[1];
sx q[1];
rz(-1.8677999) q[1];
sx q[1];
rz(-2.4261977) q[1];
rz(2.4433171) q[3];
sx q[3];
rz(-1.0812757) q[3];
sx q[3];
rz(1.5582635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58671826) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(0.47373104) q[2];
rz(2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2471152) q[0];
sx q[0];
rz(-3.101427) q[0];
sx q[0];
rz(2.1263057) q[0];
rz(2.9292551) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-2.7755348) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233519) q[0];
sx q[0];
rz(-2.2466772) q[0];
sx q[0];
rz(-0.9267207) q[0];
rz(-pi) q[1];
rz(-1.7160077) q[2];
sx q[2];
rz(-1.337707) q[2];
sx q[2];
rz(-2.9200302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3234933) q[1];
sx q[1];
rz(-2.202569) q[1];
sx q[1];
rz(-1.1786103) q[1];
x q[2];
rz(-0.29764953) q[3];
sx q[3];
rz(-0.98778703) q[3];
sx q[3];
rz(-1.7109554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4801243) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(2.4465731) q[2];
rz(-0.85092893) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(1.2581576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-2.8261321) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.618914) q[1];
sx q[1];
rz(-0.42627898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3836455) q[0];
sx q[0];
rz(-1.28048) q[0];
sx q[0];
rz(0.84602543) q[0];
rz(-2.2441909) q[2];
sx q[2];
rz(-0.36172141) q[2];
sx q[2];
rz(0.42332403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6719088) q[1];
sx q[1];
rz(-1.6571181) q[1];
sx q[1];
rz(2.0071908) q[1];
rz(-2.8180646) q[3];
sx q[3];
rz(-1.1051205) q[3];
sx q[3];
rz(2.9126008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2686501) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(-3.1128913) q[2];
rz(-0.21102333) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92006224) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(-2.9826214) q[0];
rz(-0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(1.5370625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66012525) q[0];
sx q[0];
rz(-0.9430389) q[0];
sx q[0];
rz(0.56076903) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9150429) q[2];
sx q[2];
rz(-1.2249399) q[2];
sx q[2];
rz(-2.648022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75204471) q[1];
sx q[1];
rz(-1.5524781) q[1];
sx q[1];
rz(-1.0059898) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69055478) q[3];
sx q[3];
rz(-0.46633807) q[3];
sx q[3];
rz(0.40377221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0663466) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-0.2839655) q[3];
sx q[3];
rz(0.56104463) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905106) q[0];
sx q[0];
rz(-0.53618479) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(-0.74288145) q[1];
sx q[1];
rz(-1.0123092) q[1];
sx q[1];
rz(-1.2339309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260516) q[0];
sx q[0];
rz(-1.7103467) q[0];
sx q[0];
rz(1.498879) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.874159) q[2];
sx q[2];
rz(-1.1743011) q[2];
sx q[2];
rz(-2.1934137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18542087) q[1];
sx q[1];
rz(-2.2347576) q[1];
sx q[1];
rz(-2.0645622) q[1];
rz(-pi) q[2];
rz(0.59071006) q[3];
sx q[3];
rz(-1.241467) q[3];
sx q[3];
rz(-1.2423837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10741216) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-0.91599715) q[2];
rz(-2.8902174) q[3];
sx q[3];
rz(-2.1706457) q[3];
sx q[3];
rz(2.796252) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4995572) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(0.04059759) q[0];
rz(1.9675072) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(-1.0601128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745251) q[0];
sx q[0];
rz(-2.7481348) q[0];
sx q[0];
rz(2.8889546) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7383582) q[2];
sx q[2];
rz(-3.02387) q[2];
sx q[2];
rz(0.37230834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1730949) q[1];
sx q[1];
rz(-1.8827594) q[1];
sx q[1];
rz(-0.63271823) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67146639) q[3];
sx q[3];
rz(-2.0196303) q[3];
sx q[3];
rz(-1.1374813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4982831) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(-2.3504284) q[2];
rz(-1.6061973) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(-2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20752792) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(0.15765634) q[0];
rz(-0.12570307) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(-2.8368565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3476815) q[0];
sx q[0];
rz(-1.2448989) q[0];
sx q[0];
rz(3.053717) q[0];
rz(-1.1515806) q[2];
sx q[2];
rz(-2.1682924) q[2];
sx q[2];
rz(1.0349764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3379537) q[1];
sx q[1];
rz(-1.5684761) q[1];
sx q[1];
rz(-3.0073037) q[1];
x q[2];
rz(0.79873884) q[3];
sx q[3];
rz(-1.8336356) q[3];
sx q[3];
rz(1.9376439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8704845) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(-1.452272) q[2];
rz(-2.3952386) q[3];
sx q[3];
rz(-1.690026) q[3];
sx q[3];
rz(0.11317429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8993503) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(-0.53264701) q[0];
rz(-0.38231725) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(-0.16758448) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3574651) q[0];
sx q[0];
rz(-1.7622158) q[0];
sx q[0];
rz(-3.0166059) q[0];
rz(-pi) q[1];
rz(-2.1176867) q[2];
sx q[2];
rz(-1.3452072) q[2];
sx q[2];
rz(2.0187261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51983143) q[1];
sx q[1];
rz(-1.863136) q[1];
sx q[1];
rz(-1.2755434) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83491171) q[3];
sx q[3];
rz(-1.8855699) q[3];
sx q[3];
rz(1.9251458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31183895) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(-0.46693841) q[2];
rz(-1.1030819) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5655831) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(2.3570428) q[1];
sx q[1];
rz(-0.34930925) q[1];
sx q[1];
rz(0.92225155) q[1];
rz(2.4541773) q[2];
sx q[2];
rz(-0.981642) q[2];
sx q[2];
rz(-2.8264075) q[2];
rz(1.9292694) q[3];
sx q[3];
rz(-1.8239106) q[3];
sx q[3];
rz(1.9700005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
