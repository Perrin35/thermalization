OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3990134) q[0];
sx q[0];
rz(-1.1630381) q[0];
sx q[0];
rz(-1.1385588) q[0];
rz(-0.63602716) q[1];
sx q[1];
rz(-0.50995246) q[1];
sx q[1];
rz(0.62503254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33390663) q[0];
sx q[0];
rz(-1.7117097) q[0];
sx q[0];
rz(-1.6150835) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0306084) q[2];
sx q[2];
rz(-1.5531047) q[2];
sx q[2];
rz(2.4938581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4932517) q[1];
sx q[1];
rz(-2.1025755) q[1];
sx q[1];
rz(-0.48303034) q[1];
rz(-pi) q[2];
rz(-0.81283855) q[3];
sx q[3];
rz(-1.715797) q[3];
sx q[3];
rz(-2.6036711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2193489) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(-2.7395978) q[2];
rz(0.829202) q[3];
sx q[3];
rz(-1.821937) q[3];
sx q[3];
rz(-1.7792262) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336693) q[0];
sx q[0];
rz(-0.55183721) q[0];
sx q[0];
rz(-1.1056939) q[0];
rz(0.78877527) q[1];
sx q[1];
rz(-2.5194247) q[1];
sx q[1];
rz(2.6355991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6405341) q[0];
sx q[0];
rz(-1.6468684) q[0];
sx q[0];
rz(-0.23573689) q[0];
rz(-pi) q[1];
rz(0.80319941) q[2];
sx q[2];
rz(-2.1888852) q[2];
sx q[2];
rz(1.0206285) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92365328) q[1];
sx q[1];
rz(-2.7193177) q[1];
sx q[1];
rz(2.1417066) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8333767) q[3];
sx q[3];
rz(-2.0114813) q[3];
sx q[3];
rz(-2.6868771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2371696) q[2];
sx q[2];
rz(-1.0639031) q[2];
sx q[2];
rz(-2.2696631) q[2];
rz(1.0540086) q[3];
sx q[3];
rz(-1.0796615) q[3];
sx q[3];
rz(-1.4409298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089040861) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(1.22714) q[0];
rz(0.50651208) q[1];
sx q[1];
rz(-0.7917234) q[1];
sx q[1];
rz(-0.92481771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8313356) q[0];
sx q[0];
rz(-1.4228205) q[0];
sx q[0];
rz(-0.10338115) q[0];
rz(-pi) q[1];
rz(-2.6854555) q[2];
sx q[2];
rz(-1.6326687) q[2];
sx q[2];
rz(-1.8707239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4329728) q[1];
sx q[1];
rz(-2.4390894) q[1];
sx q[1];
rz(2.0496416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1919349) q[3];
sx q[3];
rz(-0.94516813) q[3];
sx q[3];
rz(-1.715884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0792803) q[2];
sx q[2];
rz(-1.5996876) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(1.371572) q[3];
sx q[3];
rz(-2.3973231) q[3];
sx q[3];
rz(-1.3220968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1247193) q[0];
sx q[0];
rz(-1.158241) q[0];
sx q[0];
rz(0.88791263) q[0];
rz(-1.7601695) q[1];
sx q[1];
rz(-1.0180749) q[1];
sx q[1];
rz(2.6240614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84758112) q[0];
sx q[0];
rz(-1.8318307) q[0];
sx q[0];
rz(0.014495975) q[0];
rz(0.58184718) q[2];
sx q[2];
rz(-2.3334425) q[2];
sx q[2];
rz(-0.86607546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0486517) q[1];
sx q[1];
rz(-0.098628672) q[1];
sx q[1];
rz(-1.6706549) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66756396) q[3];
sx q[3];
rz(-0.57355155) q[3];
sx q[3];
rz(-1.9907601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.069328221) q[2];
sx q[2];
rz(-2.0158517) q[2];
sx q[2];
rz(2.1605055) q[2];
rz(2.6632994) q[3];
sx q[3];
rz(-0.35224733) q[3];
sx q[3];
rz(0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1961871) q[0];
sx q[0];
rz(-1.4237175) q[0];
sx q[0];
rz(0.032935306) q[0];
rz(-1.5252652) q[1];
sx q[1];
rz(-2.567629) q[1];
sx q[1];
rz(-0.93564916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0230377) q[0];
sx q[0];
rz(-1.584125) q[0];
sx q[0];
rz(2.6715735) q[0];
rz(-1.0887371) q[2];
sx q[2];
rz(-2.3080359) q[2];
sx q[2];
rz(2.478046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0914382) q[1];
sx q[1];
rz(-2.3622516) q[1];
sx q[1];
rz(0.088717566) q[1];
x q[2];
rz(2.7838681) q[3];
sx q[3];
rz(-2.5265363) q[3];
sx q[3];
rz(2.0004724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66964883) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(-1.6432537) q[2];
rz(-0.60254997) q[3];
sx q[3];
rz(-2.3157401) q[3];
sx q[3];
rz(1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5947386) q[0];
sx q[0];
rz(-2.8475519) q[0];
sx q[0];
rz(-0.038851693) q[0];
rz(-2.7815869) q[1];
sx q[1];
rz(-1.4119166) q[1];
sx q[1];
rz(0.14883277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1921944) q[0];
sx q[0];
rz(-1.4440876) q[0];
sx q[0];
rz(1.4437946) q[0];
rz(-pi) q[1];
rz(-2.932933) q[2];
sx q[2];
rz(-2.5620915) q[2];
sx q[2];
rz(2.8148871) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7578797) q[1];
sx q[1];
rz(-1.0812909) q[1];
sx q[1];
rz(2.4361324) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43977387) q[3];
sx q[3];
rz(-2.5086702) q[3];
sx q[3];
rz(-1.1732111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0845324) q[2];
sx q[2];
rz(-0.82814211) q[2];
sx q[2];
rz(-0.35489902) q[2];
rz(-2.0830294) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(2.6088349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232089) q[0];
sx q[0];
rz(-1.4223149) q[0];
sx q[0];
rz(2.3329155) q[0];
rz(-3.0478364) q[1];
sx q[1];
rz(-2.3010727) q[1];
sx q[1];
rz(0.41059986) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6855934) q[0];
sx q[0];
rz(-2.5846722) q[0];
sx q[0];
rz(2.3310175) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33538474) q[2];
sx q[2];
rz(-2.460157) q[2];
sx q[2];
rz(2.6449056) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0714662) q[1];
sx q[1];
rz(-1.1845329) q[1];
sx q[1];
rz(-1.1995308) q[1];
rz(-pi) q[2];
rz(1.4491354) q[3];
sx q[3];
rz(-0.9100737) q[3];
sx q[3];
rz(2.4796951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9984596) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(-0.03579363) q[2];
rz(-1.834747) q[3];
sx q[3];
rz(-1.1467609) q[3];
sx q[3];
rz(-0.80865639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3150113) q[0];
sx q[0];
rz(-2.1666574) q[0];
sx q[0];
rz(0.63631979) q[0];
rz(-0.66337216) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(-0.62796193) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54701824) q[0];
sx q[0];
rz(-1.5994521) q[0];
sx q[0];
rz(-3.0648607) q[0];
rz(-pi) q[1];
rz(0.3971252) q[2];
sx q[2];
rz(-1.4897122) q[2];
sx q[2];
rz(2.9604767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1455974) q[1];
sx q[1];
rz(-0.72628262) q[1];
sx q[1];
rz(2.4895345) q[1];
rz(-pi) q[2];
rz(-1.5359587) q[3];
sx q[3];
rz(-1.758134) q[3];
sx q[3];
rz(-3.0288937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58225727) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(-1.4165233) q[2];
rz(-1.0995809) q[3];
sx q[3];
rz(-1.0397592) q[3];
sx q[3];
rz(-2.288868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.002554) q[0];
sx q[0];
rz(-2.250493) q[0];
sx q[0];
rz(-2.52676) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-2.1891687) q[1];
sx q[1];
rz(2.8654548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50867453) q[0];
sx q[0];
rz(-0.44692293) q[0];
sx q[0];
rz(1.267295) q[0];
rz(-pi) q[1];
rz(1.040784) q[2];
sx q[2];
rz(-1.4333087) q[2];
sx q[2];
rz(-1.9078209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.975229) q[1];
sx q[1];
rz(-0.97492906) q[1];
sx q[1];
rz(1.0129376) q[1];
x q[2];
rz(-1.1405892) q[3];
sx q[3];
rz(-0.70616841) q[3];
sx q[3];
rz(-0.28751349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95163661) q[2];
sx q[2];
rz(-2.1026244) q[2];
sx q[2];
rz(3.0432126) q[2];
rz(-1.240587) q[3];
sx q[3];
rz(-1.0616579) q[3];
sx q[3];
rz(3.0978751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28546178) q[0];
sx q[0];
rz(-1.0482482) q[0];
sx q[0];
rz(0.35828006) q[0];
rz(-1.0048535) q[1];
sx q[1];
rz(-2.256911) q[1];
sx q[1];
rz(-2.7994432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.005363883) q[0];
sx q[0];
rz(-2.507302) q[0];
sx q[0];
rz(2.7691288) q[0];
x q[1];
rz(1.0295632) q[2];
sx q[2];
rz(-0.65766739) q[2];
sx q[2];
rz(2.5785411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6498211) q[1];
sx q[1];
rz(-0.58204466) q[1];
sx q[1];
rz(-0.19606049) q[1];
rz(-1.8588641) q[3];
sx q[3];
rz(-2.7917842) q[3];
sx q[3];
rz(1.4427078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3442012) q[2];
sx q[2];
rz(-0.48365334) q[2];
sx q[2];
rz(-0.3717711) q[2];
rz(-2.9504635) q[3];
sx q[3];
rz(-1.1647977) q[3];
sx q[3];
rz(-0.42310664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38901781) q[0];
sx q[0];
rz(-2.6448463) q[0];
sx q[0];
rz(-0.65506558) q[0];
rz(2.7816506) q[1];
sx q[1];
rz(-1.5690201) q[1];
sx q[1];
rz(1.5108861) q[1];
rz(-1.7689479) q[2];
sx q[2];
rz(-2.8336278) q[2];
sx q[2];
rz(-2.8358493) q[2];
rz(0.015040811) q[3];
sx q[3];
rz(-0.81285601) q[3];
sx q[3];
rz(2.6668919) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
