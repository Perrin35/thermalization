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
rz(2.8948957) q[0];
rz(-2.2973581) q[1];
sx q[1];
rz(-1.8480453) q[1];
sx q[1];
rz(-1.3026613) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5555252) q[0];
sx q[0];
rz(-1.925615) q[0];
sx q[0];
rz(2.8321502) q[0];
rz(-0.50452043) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(-0.96790403) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0595536) q[1];
sx q[1];
rz(-0.63602018) q[1];
sx q[1];
rz(0.18853699) q[1];
x q[2];
rz(2.2095895) q[3];
sx q[3];
rz(-2.9607256) q[3];
sx q[3];
rz(2.7993134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.3641597) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(1.9470661) q[2];
rz(-1.2398237) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(-1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251004) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-2.4745353) q[0];
rz(-2.8705813) q[1];
sx q[1];
rz(-1.4457116) q[1];
sx q[1];
rz(-2.0858696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0843029) q[0];
sx q[0];
rz(-2.192756) q[0];
sx q[0];
rz(-0.33226407) q[0];
rz(-pi) q[1];
rz(1.5178568) q[2];
sx q[2];
rz(-0.38489562) q[2];
sx q[2];
rz(3.0786849) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70433863) q[1];
sx q[1];
rz(-1.9074215) q[1];
sx q[1];
rz(-3.047154) q[1];
rz(2.3930644) q[3];
sx q[3];
rz(-1.0077884) q[3];
sx q[3];
rz(-2.9070118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6229652) q[2];
sx q[2];
rz(-2.3845606) q[2];
sx q[2];
rz(0.61678994) q[2];
rz(-2.8957497) q[3];
sx q[3];
rz(-1.5522233) q[3];
sx q[3];
rz(1.2639379) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(-0.47897419) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(-1.3165547) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265357) q[0];
sx q[0];
rz(-2.2000072) q[0];
sx q[0];
rz(0.75124426) q[0];
rz(-pi) q[1];
rz(1.1321849) q[2];
sx q[2];
rz(-0.45044611) q[2];
sx q[2];
rz(0.23368719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41539792) q[1];
sx q[1];
rz(-1.2737927) q[1];
sx q[1];
rz(-0.71539496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69200055) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58671826) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(-2.6678616) q[2];
rz(-2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(1.0152869) q[0];
rz(-2.9292551) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(2.7755348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.191712) q[0];
sx q[0];
rz(-2.0583344) q[0];
sx q[0];
rz(-0.78678188) q[0];
x q[1];
rz(0.54738657) q[2];
sx q[2];
rz(-0.27392188) q[2];
sx q[2];
rz(0.78597921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8541222) q[1];
sx q[1];
rz(-0.72919289) q[1];
sx q[1];
rz(-2.6602938) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9892392) q[3];
sx q[3];
rz(-0.64665808) q[3];
sx q[3];
rz(1.2025976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(0.69501957) q[2];
rz(-2.2906637) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.2581576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51915351) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-2.8261321) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-2.7153137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56458679) q[0];
sx q[0];
rz(-2.2590911) q[0];
sx q[0];
rz(-2.7618963) q[0];
rz(2.909864) q[2];
sx q[2];
rz(-1.2905057) q[2];
sx q[2];
rz(2.0120399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91831827) q[1];
sx q[1];
rz(-0.4443109) q[1];
sx q[1];
rz(1.3688509) q[1];
rz(-pi) q[2];
rz(-1.0067389) q[3];
sx q[3];
rz(-2.5814179) q[3];
sx q[3];
rz(2.2711846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2686501) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(-0.028701393) q[2];
rz(-0.21102333) q[3];
sx q[3];
rz(-0.70091453) q[3];
sx q[3];
rz(2.4215462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(0.15897121) q[0];
rz(2.4863415) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(1.5370625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814674) q[0];
sx q[0];
rz(-0.9430389) q[0];
sx q[0];
rz(0.56076903) q[0];
x q[1];
rz(0.36559029) q[2];
sx q[2];
rz(-1.2477008) q[2];
sx q[2];
rz(-1.9434203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3517396) q[1];
sx q[1];
rz(-2.5765214) q[1];
sx q[1];
rz(1.6050102) q[1];
rz(-pi) q[2];
rz(1.8810684) q[3];
sx q[3];
rz(-1.9247484) q[3];
sx q[3];
rz(-1.1502532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0752461) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(-3.1250478) q[3];
sx q[3];
rz(-0.2839655) q[3];
sx q[3];
rz(2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(-0.91019994) q[0];
rz(2.3987112) q[1];
sx q[1];
rz(-1.0123092) q[1];
sx q[1];
rz(1.9076617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5963579) q[0];
sx q[0];
rz(-1.6420134) q[0];
sx q[0];
rz(-0.13990732) q[0];
x q[1];
rz(-0.26743369) q[2];
sx q[2];
rz(-1.9672915) q[2];
sx q[2];
rz(-2.1934137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2382395) q[1];
sx q[1];
rz(-2.337114) q[1];
sx q[1];
rz(-2.5969905) q[1];
rz(-pi) q[2];
rz(0.59071006) q[3];
sx q[3];
rz(-1.241467) q[3];
sx q[3];
rz(1.8992089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-0.91599715) q[2];
rz(0.25137526) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995572) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-3.1009951) q[0];
rz(1.9675072) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(2.0814799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6396234) q[0];
sx q[0];
rz(-1.9511128) q[0];
sx q[0];
rz(-1.4674076) q[0];
x q[1];
rz(-1.6868851) q[2];
sx q[2];
rz(-1.5512067) q[2];
sx q[2];
rz(-1.7766812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6190784) q[1];
sx q[1];
rz(-0.97303094) q[1];
sx q[1];
rz(-1.1903711) q[1];
x q[2];
rz(2.4701263) q[3];
sx q[3];
rz(-1.1219624) q[3];
sx q[3];
rz(2.0041114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64330953) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(2.3504284) q[2];
rz(1.5353954) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(-0.15765634) q[0];
rz(0.12570307) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(-2.8368565) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9466772) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(-1.8978682) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9900121) q[2];
sx q[2];
rz(-2.1682924) q[2];
sx q[2];
rz(-2.1066163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.803639) q[1];
sx q[1];
rz(-1.5731166) q[1];
sx q[1];
rz(0.13428899) q[1];
rz(2.7823506) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(-0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2711082) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.452272) q[2];
rz(2.3952386) q[3];
sx q[3];
rz(-1.690026) q[3];
sx q[3];
rz(3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8993503) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(-0.53264701) q[0];
rz(-2.7592754) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(0.16758448) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9411128) q[0];
sx q[0];
rz(-0.22819209) q[0];
sx q[0];
rz(-2.1424293) q[0];
rz(-1.0239059) q[2];
sx q[2];
rz(-1.7963855) q[2];
sx q[2];
rz(2.0187261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1780562) q[1];
sx q[1];
rz(-1.85317) q[1];
sx q[1];
rz(-2.8368188) q[1];
rz(-pi) q[2];
rz(-2.7277023) q[3];
sx q[3];
rz(-0.87855708) q[3];
sx q[3];
rz(0.62769753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(0.46693841) q[2];
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
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5655831) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(0.78454984) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(-2.2837737) q[2];
sx q[2];
rz(-2.1265278) q[2];
sx q[2];
rz(2.3139755) q[2];
rz(0.26950027) q[3];
sx q[3];
rz(-1.2242347) q[3];
sx q[3];
rz(0.4927529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
