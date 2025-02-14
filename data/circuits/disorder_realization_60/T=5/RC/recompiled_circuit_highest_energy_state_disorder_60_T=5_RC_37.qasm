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
rz(-2.2973581) q[1];
sx q[1];
rz(-1.8480453) q[1];
sx q[1];
rz(1.8389314) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58606746) q[0];
sx q[0];
rz(-1.925615) q[0];
sx q[0];
rz(-0.3094425) q[0];
rz(-2.6370722) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(-2.1736886) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.0820391) q[1];
sx q[1];
rz(-0.63602018) q[1];
sx q[1];
rz(-0.18853699) q[1];
x q[2];
rz(-2.2095895) q[3];
sx q[3];
rz(-0.18086704) q[3];
sx q[3];
rz(2.7993134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(-1.9470661) q[2];
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
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164923) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(0.66705739) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.4457116) q[1];
sx q[1];
rz(-1.0557231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843029) q[0];
sx q[0];
rz(-2.192756) q[0];
sx q[0];
rz(-2.8093286) q[0];
rz(-3.1201601) q[2];
sx q[2];
rz(-1.9551245) q[2];
sx q[2];
rz(-3.1357946) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70433863) q[1];
sx q[1];
rz(-1.9074215) q[1];
sx q[1];
rz(-3.047154) q[1];
rz(-pi) q[2];
rz(0.74774489) q[3];
sx q[3];
rz(-0.90255957) q[3];
sx q[3];
rz(-1.2838319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51862741) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(2.5248027) q[2];
rz(-0.24584298) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.2639379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-2.9611294) q[0];
rz(0.47897419) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(1.825038) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1235629) q[0];
sx q[0];
rz(-0.93864894) q[0];
sx q[0];
rz(-2.324047) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9834824) q[2];
sx q[2];
rz(-1.7567593) q[2];
sx q[2];
rz(-0.93753147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6613343) q[1];
sx q[1];
rz(-0.76441718) q[1];
sx q[1];
rz(-2.7050128) q[1];
rz(-0.69827551) q[3];
sx q[3];
rz(-1.0812757) q[3];
sx q[3];
rz(-1.5833291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5548744) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(-2.6678616) q[2];
rz(2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.1562851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(-1.0152869) q[0];
rz(0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-0.36605787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9498807) q[0];
sx q[0];
rz(-1.0832583) q[0];
sx q[0];
rz(-0.78678188) q[0];
rz(1.7160077) q[2];
sx q[2];
rz(-1.8038857) q[2];
sx q[2];
rz(0.22156246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4868813) q[1];
sx q[1];
rz(-1.8843448) q[1];
sx q[1];
rz(-2.4717719) q[1];
rz(-pi) q[2];
rz(-1.1523535) q[3];
sx q[3];
rz(-2.4949346) q[3];
sx q[3];
rz(1.2025976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-2.4465731) q[2];
rz(-2.2906637) q[3];
sx q[3];
rz(-1.3635819) q[3];
sx q[3];
rz(-1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-1.618914) q[1];
sx q[1];
rz(-0.42627898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75794719) q[0];
sx q[0];
rz(-1.8611127) q[0];
sx q[0];
rz(2.2955672) q[0];
rz(-pi) q[1];
rz(-1.8583723) q[2];
sx q[2];
rz(-1.3482665) q[2];
sx q[2];
rz(2.6351647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46968383) q[1];
sx q[1];
rz(-1.6571181) q[1];
sx q[1];
rz(2.0071908) q[1];
rz(-pi) q[2];
rz(1.0834094) q[3];
sx q[3];
rz(-1.8588239) q[3];
sx q[3];
rz(1.9492287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8729426) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(-0.028701393) q[2];
rz(-2.9305693) q[3];
sx q[3];
rz(-0.70091453) q[3];
sx q[3];
rz(0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92006224) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(2.9826214) q[0];
rz(-2.4863415) q[1];
sx q[1];
rz(-2.7924004) q[1];
sx q[1];
rz(1.5370625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814674) q[0];
sx q[0];
rz(-2.1985538) q[0];
sx q[0];
rz(0.56076903) q[0];
rz(-pi) q[1];
rz(-2.3889306) q[2];
sx q[2];
rz(-0.48303451) q[2];
sx q[2];
rz(0.31980425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75204471) q[1];
sx q[1];
rz(-1.5524781) q[1];
sx q[1];
rz(-1.0059898) q[1];
rz(-0.37015583) q[3];
sx q[3];
rz(-1.8612544) q[3];
sx q[3];
rz(-2.6103717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0663466) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(0.1417024) q[2];
rz(-3.1250478) q[3];
sx q[3];
rz(-0.2839655) q[3];
sx q[3];
rz(2.580548) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510821) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(-2.2313927) q[0];
rz(0.74288145) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(1.9076617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260516) q[0];
sx q[0];
rz(-1.7103467) q[0];
sx q[0];
rz(-1.6427137) q[0];
x q[1];
rz(-0.26743369) q[2];
sx q[2];
rz(-1.9672915) q[2];
sx q[2];
rz(-2.1934137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9561718) q[1];
sx q[1];
rz(-2.2347576) q[1];
sx q[1];
rz(-2.0645622) q[1];
rz(-2.5508826) q[3];
sx q[3];
rz(-1.241467) q[3];
sx q[3];
rz(1.8992089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(2.2255955) q[2];
rz(0.25137526) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6420355) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-0.04059759) q[0];
rz(1.9675072) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(-2.0814799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019693) q[0];
sx q[0];
rz(-1.1904799) q[0];
sx q[0];
rz(1.4674076) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7383582) q[2];
sx q[2];
rz(-0.11772269) q[2];
sx q[2];
rz(2.7692843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9684978) q[1];
sx q[1];
rz(-1.8827594) q[1];
sx q[1];
rz(-0.63271823) q[1];
rz(-1.0193018) q[3];
sx q[3];
rz(-2.1658033) q[3];
sx q[3];
rz(2.3762356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64330953) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(-2.3504284) q[2];
rz(-1.5353954) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(0.78996381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.9340647) q[0];
sx q[0];
rz(-2.9495033) q[0];
sx q[0];
rz(-0.15765634) q[0];
rz(0.12570307) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(-2.8368565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5253741) q[0];
sx q[0];
rz(-2.8044639) q[0];
sx q[0];
rz(1.824877) q[0];
rz(-0.64027787) q[2];
sx q[2];
rz(-1.9140179) q[2];
sx q[2];
rz(2.8514112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.803639) q[1];
sx q[1];
rz(-1.5684761) q[1];
sx q[1];
rz(-3.0073037) q[1];
x q[2];
rz(-2.7823506) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(-3.0225282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8704845) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.452272) q[2];
rz(0.74635402) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-0.11317429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24224237) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(0.53264701) q[0];
rz(-0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(0.16758448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2004799) q[0];
sx q[0];
rz(-0.22819209) q[0];
sx q[0];
rz(2.1424293) q[0];
rz(1.1551935) q[2];
sx q[2];
rz(-2.5544082) q[2];
sx q[2];
rz(-0.09584643) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51983143) q[1];
sx q[1];
rz(-1.863136) q[1];
sx q[1];
rz(1.8660493) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7277023) q[3];
sx q[3];
rz(-2.2630356) q[3];
sx q[3];
rz(-2.5138951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(0.46693841) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.1934049) q[3];
sx q[3];
rz(-1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5655831) q[0];
sx q[0];
rz(-0.82721114) q[0];
sx q[0];
rz(0.4785434) q[0];
rz(0.78454984) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(-0.68741531) q[2];
sx q[2];
rz(-0.981642) q[2];
sx q[2];
rz(-2.8264075) q[2];
rz(-1.9292694) q[3];
sx q[3];
rz(-1.3176821) q[3];
sx q[3];
rz(-1.1715922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
