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
rz(-1.3026613) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3299574) q[0];
sx q[0];
rz(-0.46644743) q[0];
sx q[0];
rz(0.88282013) q[0];
rz(-pi) q[1];
rz(2.6370722) q[2];
sx q[2];
rz(-2.2103643) q[2];
sx q[2];
rz(0.96790403) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0820391) q[1];
sx q[1];
rz(-0.63602018) q[1];
sx q[1];
rz(2.9530557) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2095895) q[3];
sx q[3];
rz(-0.18086704) q[3];
sx q[3];
rz(0.34227926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(-1.1945266) q[2];
rz(1.2398237) q[3];
sx q[3];
rz(-1.4907962) q[3];
sx q[3];
rz(-1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251004) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-2.4745353) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0843029) q[0];
sx q[0];
rz(-2.192756) q[0];
sx q[0];
rz(2.8093286) q[0];
rz(-pi) q[1];
rz(1.5178568) q[2];
sx q[2];
rz(-0.38489562) q[2];
sx q[2];
rz(-0.062907779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.437254) q[1];
sx q[1];
rz(-1.2341712) q[1];
sx q[1];
rz(3.047154) q[1];
rz(-pi) q[2];
rz(-0.74774489) q[3];
sx q[3];
rz(-0.90255957) q[3];
sx q[3];
rz(1.2838319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51862741) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(0.24584298) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-0.73391947) q[0];
sx q[0];
rz(2.9611294) q[0];
rz(0.47897419) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(-1.825038) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265357) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(0.75124426) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20256217) q[2];
sx q[2];
rz(-1.975946) q[2];
sx q[2];
rz(-0.7140401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2352202) q[1];
sx q[1];
rz(-2.2487469) q[1];
sx q[1];
rz(1.9560019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4433171) q[3];
sx q[3];
rz(-1.0812757) q[3];
sx q[3];
rz(1.5582635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5548744) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(0.47373104) q[2];
rz(2.6321865) q[3];
sx q[3];
rz(-1.3207057) q[3];
sx q[3];
rz(1.1562851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(2.1263057) q[0];
rz(-0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-2.7755348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057945874) q[0];
sx q[0];
rz(-0.89712954) q[0];
sx q[0];
rz(-2.4988453) q[0];
rz(1.7160077) q[2];
sx q[2];
rz(-1.8038857) q[2];
sx q[2];
rz(0.22156246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4868813) q[1];
sx q[1];
rz(-1.8843448) q[1];
sx q[1];
rz(2.4717719) q[1];
rz(-pi) q[2];
rz(-1.9892392) q[3];
sx q[3];
rz(-2.4949346) q[3];
sx q[3];
rz(1.9389951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-0.69501957) q[2];
rz(-0.85092893) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.618914) q[1];
sx q[1];
rz(2.7153137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.016092) q[0];
sx q[0];
rz(-2.3707485) q[0];
sx q[0];
rz(-1.9941814) q[0];
x q[1];
rz(2.909864) q[2];
sx q[2];
rz(-1.2905057) q[2];
sx q[2];
rz(-1.1295527) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6719088) q[1];
sx q[1];
rz(-1.4844746) q[1];
sx q[1];
rz(-2.0071908) q[1];
rz(-pi) q[2];
rz(-1.0067389) q[3];
sx q[3];
rz(-2.5814179) q[3];
sx q[3];
rz(-0.87040802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2686501) q[2];
sx q[2];
rz(-0.23317569) q[2];
sx q[2];
rz(0.028701393) q[2];
rz(2.9305693) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92006224) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(2.9826214) q[0];
rz(-0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(-1.6045301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5842954) q[0];
sx q[0];
rz(-1.125815) q[0];
sx q[0];
rz(-0.86229445) q[0];
x q[1];
rz(1.2265497) q[2];
sx q[2];
rz(-1.2249399) q[2];
sx q[2];
rz(2.648022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83035806) q[1];
sx q[1];
rz(-1.0060961) q[1];
sx q[1];
rz(0.021685251) q[1];
rz(2.4510379) q[3];
sx q[3];
rz(-0.46633807) q[3];
sx q[3];
rz(0.40377221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0752461) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(-2.3987112) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(-1.2339309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.493422) q[0];
sx q[0];
rz(-0.1568846) q[0];
sx q[0];
rz(-0.4728526) q[0];
rz(-pi) q[1];
rz(-1.1612438) q[2];
sx q[2];
rz(-1.8170333) q[2];
sx q[2];
rz(-0.72803942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90335315) q[1];
sx q[1];
rz(-2.337114) q[1];
sx q[1];
rz(-0.54460214) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5911932) q[3];
sx q[3];
rz(-2.4749651) q[3];
sx q[3];
rz(0.77778274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-2.2255955) q[2];
rz(0.25137526) q[3];
sx q[3];
rz(-2.1706457) q[3];
sx q[3];
rz(2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6420355) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-3.1009951) q[0];
rz(-1.9675072) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(-1.0601128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019693) q[0];
sx q[0];
rz(-1.1904799) q[0];
sx q[0];
rz(1.4674076) q[0];
rz(1.7383582) q[2];
sx q[2];
rz(-0.11772269) q[2];
sx q[2];
rz(0.37230834) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1730949) q[1];
sx q[1];
rz(-1.2588333) q[1];
sx q[1];
rz(-0.63271823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67146639) q[3];
sx q[3];
rz(-1.1219624) q[3];
sx q[3];
rz(-1.1374813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20752792) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(0.15765634) q[0];
rz(-3.0158896) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(0.30473614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9466772) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(1.8978682) q[0];
rz(-pi) q[1];
rz(-2.5013148) q[2];
sx q[2];
rz(-1.9140179) q[2];
sx q[2];
rz(-2.8514112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.803639) q[1];
sx q[1];
rz(-1.5731166) q[1];
sx q[1];
rz(-3.0073037) q[1];
x q[2];
rz(-2.7823506) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2711082) q[2];
sx q[2];
rz(-1.5616337) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8993503) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(2.6089456) q[0];
rz(-2.7592754) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(-2.9740082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3310247) q[0];
sx q[0];
rz(-1.448104) q[0];
sx q[0];
rz(1.3779089) q[0];
x q[1];
rz(-1.0239059) q[2];
sx q[2];
rz(-1.7963855) q[2];
sx q[2];
rz(2.0187261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1780562) q[1];
sx q[1];
rz(-1.2884226) q[1];
sx q[1];
rz(0.30477384) q[1];
rz(-pi) q[2];
rz(-0.83491171) q[3];
sx q[3];
rz(-1.8855699) q[3];
sx q[3];
rz(-1.2164468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(0.46693841) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(-1.9095437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(2.3570428) q[1];
sx q[1];
rz(-0.34930925) q[1];
sx q[1];
rz(0.92225155) q[1];
rz(-0.81132728) q[2];
sx q[2];
rz(-0.87292508) q[2];
sx q[2];
rz(1.2909918) q[2];
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
