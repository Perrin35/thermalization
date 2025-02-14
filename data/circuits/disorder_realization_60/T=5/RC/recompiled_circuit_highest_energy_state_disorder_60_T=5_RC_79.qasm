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
rz(3.9858272) q[1];
sx q[1];
rz(1.2935473) q[1];
sx q[1];
rz(10.727439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87410754) q[0];
sx q[0];
rz(-1.2812072) q[0];
sx q[0];
rz(-1.1998313) q[0];
x q[1];
rz(0.86645007) q[2];
sx q[2];
rz(-1.9690919) q[2];
sx q[2];
rz(0.9212538) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9907551) q[1];
sx q[1];
rz(-2.1937943) q[1];
sx q[1];
rz(1.7083108) q[1];
rz(-pi) q[2];
rz(-1.7165623) q[3];
sx q[3];
rz(-1.4633388) q[3];
sx q[3];
rz(-1.2821357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(-1.9470661) q[2];
rz(-1.2398237) q[3];
sx q[3];
rz(-1.4907962) q[3];
sx q[3];
rz(-1.4136081) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251004) q[0];
sx q[0];
rz(-1.1273552) q[0];
sx q[0];
rz(-0.66705739) q[0];
rz(-2.8705813) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(2.0858696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572898) q[0];
sx q[0];
rz(-2.192756) q[0];
sx q[0];
rz(-0.33226407) q[0];
rz(-pi) q[1];
rz(1.6237359) q[2];
sx q[2];
rz(-0.38489562) q[2];
sx q[2];
rz(-3.0786849) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2438587) q[1];
sx q[1];
rz(-1.6599201) q[1];
sx q[1];
rz(-1.9088163) q[1];
x q[2];
rz(0.74852826) q[3];
sx q[3];
rz(-1.0077884) q[3];
sx q[3];
rz(-0.23458086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51862741) q[2];
sx q[2];
rz(-2.3845606) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(2.8957497) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.2639379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-2.9611294) q[0];
rz(2.6626185) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(1.3165547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.115057) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(2.3903484) q[0];
rz(-1.9834824) q[2];
sx q[2];
rz(-1.7567593) q[2];
sx q[2];
rz(0.93753147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4802584) q[1];
sx q[1];
rz(-2.3771755) q[1];
sx q[1];
rz(0.43657984) q[1];
rz(-pi) q[2];
rz(0.96305029) q[3];
sx q[3];
rz(-0.96745771) q[3];
sx q[3];
rz(-0.36336366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5548744) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(2.6678616) q[2];
rz(0.50940618) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2471152) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(2.1263057) q[0];
rz(0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(2.7755348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.191712) q[0];
sx q[0];
rz(-1.0832583) q[0];
sx q[0];
rz(-0.78678188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7160077) q[2];
sx q[2];
rz(-1.8038857) q[2];
sx q[2];
rz(-0.22156246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3234933) q[1];
sx q[1];
rz(-0.93902367) q[1];
sx q[1];
rz(-1.1786103) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96694209) q[3];
sx q[3];
rz(-1.3234541) q[3];
sx q[3];
rz(3.1144457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4801243) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(-0.69501957) q[2];
rz(-2.2906637) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.2581576) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(0.31546053) q[0];
rz(-2.8568965) q[1];
sx q[1];
rz(-1.618914) q[1];
sx q[1];
rz(0.42627898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3836455) q[0];
sx q[0];
rz(-1.28048) q[0];
sx q[0];
rz(-0.84602543) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89740173) q[2];
sx q[2];
rz(-0.36172141) q[2];
sx q[2];
rz(-2.7182686) q[2];
rz(-pi) q[3];
x q[3];
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
x q[2];
rz(-2.1348537) q[3];
sx q[3];
rz(-0.56017471) q[3];
sx q[3];
rz(2.2711846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8729426) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(-0.028701393) q[2];
rz(2.9305693) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(0.72004643) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(-2.9826214) q[0];
rz(-0.65525118) q[1];
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
rz(0.55729729) q[0];
sx q[0];
rz(-1.125815) q[0];
sx q[0];
rz(-2.2792982) q[0];
rz(1.2265497) q[2];
sx q[2];
rz(-1.9166528) q[2];
sx q[2];
rz(0.49357061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78985301) q[1];
sx q[1];
rz(-0.56507128) q[1];
sx q[1];
rz(1.5365824) q[1];
x q[2];
rz(-1.8810684) q[3];
sx q[3];
rz(-1.2168443) q[3];
sx q[3];
rz(1.9913395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0752461) q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(-0.91019994) q[0];
rz(-0.74288145) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(-1.9076617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260516) q[0];
sx q[0];
rz(-1.7103467) q[0];
sx q[0];
rz(-1.6427137) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1612438) q[2];
sx q[2];
rz(-1.8170333) q[2];
sx q[2];
rz(2.4135532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0651111) q[1];
sx q[1];
rz(-1.1882768) q[1];
sx q[1];
rz(0.72648813) q[1];
x q[2];
rz(-2.5508826) q[3];
sx q[3];
rz(-1.241467) q[3];
sx q[3];
rz(-1.2423837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.567652) q[2];
sx q[2];
rz(2.2255955) q[2];
rz(-2.8902174) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(-2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995572) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(3.1009951) q[0];
rz(1.1740855) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(2.0814799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745251) q[0];
sx q[0];
rz(-0.39345783) q[0];
sx q[0];
rz(2.8889546) q[0];
rz(-pi) q[1];
rz(-1.4547075) q[2];
sx q[2];
rz(-1.5512067) q[2];
sx q[2];
rz(-1.3649114) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0012944) q[1];
sx q[1];
rz(-2.4457275) q[1];
sx q[1];
rz(-2.6423179) q[1];
rz(-pi) q[2];
rz(-0.67146639) q[3];
sx q[3];
rz(-1.1219624) q[3];
sx q[3];
rz(2.0041114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4982831) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(-0.79116428) q[2];
rz(-1.6061973) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(-0.78996381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(-2.9839363) q[0];
rz(-3.0158896) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(2.8368565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9466772) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(-1.8978682) q[0];
rz(-1.1515806) q[2];
sx q[2];
rz(-0.97330026) q[2];
sx q[2];
rz(-1.0349764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25001502) q[1];
sx q[1];
rz(-0.13430891) q[1];
sx q[1];
rz(-3.1242643) q[1];
rz(-pi) q[2];
rz(-0.79873884) q[3];
sx q[3];
rz(-1.307957) q[3];
sx q[3];
rz(-1.2039487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8704845) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.452272) q[2];
rz(-2.3952386) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-0.11317429) q[3];
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
rz(2.8993503) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(-2.6089456) q[0];
rz(2.7592754) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(-0.16758448) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3310247) q[0];
sx q[0];
rz(-1.6934886) q[0];
sx q[0];
rz(-1.3779089) q[0];
x q[1];
rz(-1.9863992) q[2];
sx q[2];
rz(-0.58718449) q[2];
sx q[2];
rz(0.09584643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51983143) q[1];
sx q[1];
rz(-1.2784567) q[1];
sx q[1];
rz(-1.2755434) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1191649) q[3];
sx q[3];
rz(-0.78863025) q[3];
sx q[3];
rz(-0.024922289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31183895) q[2];
sx q[2];
rz(-1.6646174) q[2];
sx q[2];
rz(-0.46693841) q[2];
rz(-2.0385108) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(-1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(0.78454984) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(0.81132728) q[2];
sx q[2];
rz(-2.2686676) q[2];
sx q[2];
rz(-1.8506008) q[2];
rz(-2.2060903) q[3];
sx q[3];
rz(-0.43564921) q[3];
sx q[3];
rz(2.9516006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
