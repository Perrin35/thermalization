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
rz(0.84423455) q[1];
sx q[1];
rz(-1.2935473) q[1];
sx q[1];
rz(-1.8389314) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5555252) q[0];
sx q[0];
rz(-1.925615) q[0];
sx q[0];
rz(-2.8321502) q[0];
rz(2.1470492) q[2];
sx q[2];
rz(-2.3495397) q[2];
sx q[2];
rz(0.22135569) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6410632) q[1];
sx q[1];
rz(-1.682356) q[1];
sx q[1];
rz(-0.62749933) q[1];
rz(0.93200316) q[3];
sx q[3];
rz(-2.9607256) q[3];
sx q[3];
rz(-2.7993134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.777433) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(-1.1945266) q[2];
rz(-1.901769) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(-1.4136081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251004) q[0];
sx q[0];
rz(-1.1273552) q[0];
sx q[0];
rz(-2.4745353) q[0];
rz(2.8705813) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5496084) q[0];
sx q[0];
rz(-2.4469564) q[0];
sx q[0];
rz(-1.1437835) q[0];
rz(-pi) q[1];
rz(3.1201601) q[2];
sx q[2];
rz(-1.1864682) q[2];
sx q[2];
rz(0.0057980428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7165274) q[1];
sx q[1];
rz(-2.7924573) q[1];
sx q[1];
rz(-1.3075816) q[1];
rz(-pi) q[2];
rz(0.74852826) q[3];
sx q[3];
rz(-1.0077884) q[3];
sx q[3];
rz(-0.23458086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51862741) q[2];
sx q[2];
rz(-2.3845606) q[2];
sx q[2];
rz(-2.5248027) q[2];
rz(-0.24584298) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.2639379) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-2.9611294) q[0];
rz(-0.47897419) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(1.3165547) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.115057) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(-0.75124426) q[0];
rz(-1.1321849) q[2];
sx q[2];
rz(-2.6911465) q[2];
sx q[2];
rz(0.23368719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6613343) q[1];
sx q[1];
rz(-2.3771755) q[1];
sx q[1];
rz(-0.43657984) q[1];
rz(-0.96305029) q[3];
sx q[3];
rz(-0.96745771) q[3];
sx q[3];
rz(0.36336366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5548744) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(0.47373104) q[2];
rz(2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.1562851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2471152) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(1.0152869) q[0];
rz(-0.21233755) q[1];
sx q[1];
rz(-1.5539853) q[1];
sx q[1];
rz(-0.36605787) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057945874) q[0];
sx q[0];
rz(-0.89712954) q[0];
sx q[0];
rz(-2.4988453) q[0];
x q[1];
rz(-1.7160077) q[2];
sx q[2];
rz(-1.8038857) q[2];
sx q[2];
rz(-0.22156246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3234933) q[1];
sx q[1];
rz(-2.202569) q[1];
sx q[1];
rz(1.1786103) q[1];
x q[2];
rz(-0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(-0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(0.69501957) q[2];
rz(0.85092893) q[3];
sx q[3];
rz(-1.3635819) q[3];
sx q[3];
rz(-1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-0.31546053) q[0];
rz(-2.8568965) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-0.42627898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3836455) q[0];
sx q[0];
rz(-1.8611127) q[0];
sx q[0];
rz(-2.2955672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8583723) q[2];
sx q[2];
rz(-1.3482665) q[2];
sx q[2];
rz(-2.6351647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0002928) q[1];
sx q[1];
rz(-2.0054549) q[1];
sx q[1];
rz(3.0463957) q[1];
rz(-pi) q[2];
rz(1.0834094) q[3];
sx q[3];
rz(-1.2827687) q[3];
sx q[3];
rz(-1.9492287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(-0.15897121) q[0];
rz(2.4863415) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(-1.6045301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55729729) q[0];
sx q[0];
rz(-2.0157776) q[0];
sx q[0];
rz(0.86229445) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2265497) q[2];
sx q[2];
rz(-1.2249399) q[2];
sx q[2];
rz(-0.49357061) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.4510379) q[3];
sx q[3];
rz(-0.46633807) q[3];
sx q[3];
rz(0.40377221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0752461) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(0.1417024) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(-0.56104463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510821) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(-2.3987112) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(1.9076617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6481706) q[0];
sx q[0];
rz(-2.9847081) q[0];
sx q[0];
rz(-2.6687401) q[0];
rz(-pi) q[1];
rz(-2.874159) q[2];
sx q[2];
rz(-1.1743011) q[2];
sx q[2];
rz(0.948179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0764816) q[1];
sx q[1];
rz(-1.9533159) q[1];
sx q[1];
rz(-0.72648813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5508826) q[3];
sx q[3];
rz(-1.9001257) q[3];
sx q[3];
rz(-1.8992089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10741216) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(2.2255955) q[2];
rz(2.8902174) q[3];
sx q[3];
rz(-2.1706457) q[3];
sx q[3];
rz(0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6420355) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-0.04059759) q[0];
rz(1.1740855) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(1.0601128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1112633) q[0];
sx q[0];
rz(-1.4748186) q[0];
sx q[0];
rz(-2.7594271) q[0];
rz(-3.1218704) q[2];
sx q[2];
rz(-1.6868627) q[2];
sx q[2];
rz(-2.937992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1402982) q[1];
sx q[1];
rz(-2.4457275) q[1];
sx q[1];
rz(-0.49927478) q[1];
rz(0.65877093) q[3];
sx q[3];
rz(-0.78785721) q[3];
sx q[3];
rz(0.066490563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(2.3504284) q[2];
rz(-1.5353954) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(-0.15765634) q[0];
rz(3.0158896) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(0.30473614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7939111) q[0];
sx q[0];
rz(-1.8966938) q[0];
sx q[0];
rz(3.053717) q[0];
x q[1];
rz(-0.53908277) q[2];
sx q[2];
rz(-0.7149019) q[2];
sx q[2];
rz(-1.7049005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.803639) q[1];
sx q[1];
rz(-1.5731166) q[1];
sx q[1];
rz(3.0073037) q[1];
rz(0.79873884) q[3];
sx q[3];
rz(-1.307957) q[3];
sx q[3];
rz(1.2039487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8704845) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(-1.6893207) q[2];
rz(2.3952386) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24224237) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(-0.53264701) q[0];
rz(0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(-0.16758448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2004799) q[0];
sx q[0];
rz(-0.22819209) q[0];
sx q[0];
rz(2.1424293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8791076) q[2];
sx q[2];
rz(-1.0392611) q[2];
sx q[2];
rz(-2.558311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96353647) q[1];
sx q[1];
rz(-1.85317) q[1];
sx q[1];
rz(0.30477384) q[1];
rz(1.1191649) q[3];
sx q[3];
rz(-2.3529624) q[3];
sx q[3];
rz(-0.024922289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.3302654) q[2];
sx q[2];
rz(-2.2686676) q[2];
sx q[2];
rz(-1.8506008) q[2];
rz(1.2123232) q[3];
sx q[3];
rz(-1.3176821) q[3];
sx q[3];
rz(-1.1715922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
