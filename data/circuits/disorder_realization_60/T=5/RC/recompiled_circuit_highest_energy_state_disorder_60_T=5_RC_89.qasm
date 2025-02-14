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
rz(1.4580392) q[0];
sx q[0];
rz(-0.91705051) q[0];
sx q[0];
rz(0.24669692) q[0];
rz(-2.2973581) q[1];
sx q[1];
rz(-1.8480453) q[1];
sx q[1];
rz(-1.3026613) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674851) q[0];
sx q[0];
rz(-1.2812072) q[0];
sx q[0];
rz(-1.1998313) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6370722) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(0.96790403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0820391) q[1];
sx q[1];
rz(-2.5055725) q[1];
sx q[1];
rz(2.9530557) q[1];
rz(-pi) q[2];
rz(3.0329923) q[3];
sx q[3];
rz(-1.7157156) q[3];
sx q[3];
rz(0.30440457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.777433) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(1.1945266) q[2];
rz(1.901769) q[3];
sx q[3];
rz(-1.4907962) q[3];
sx q[3];
rz(1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251004) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(2.4745353) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5496084) q[0];
sx q[0];
rz(-0.69463629) q[0];
sx q[0];
rz(1.9978092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9552043) q[2];
sx q[2];
rz(-1.5509275) q[2];
sx q[2];
rz(1.5846313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7165274) q[1];
sx q[1];
rz(-0.34913539) q[1];
sx q[1];
rz(1.8340111) q[1];
rz(-pi) q[2];
rz(2.3930644) q[3];
sx q[3];
rz(-2.1338042) q[3];
sx q[3];
rz(-0.23458086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701876) q[0];
sx q[0];
rz(-0.73391947) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(2.6626185) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(1.825038) q[1];
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
rz(0.75124426) q[0];
rz(-pi) q[1];
rz(1.1321849) q[2];
sx q[2];
rz(-2.6911465) q[2];
sx q[2];
rz(2.9079055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7261947) q[1];
sx q[1];
rz(-1.2737927) q[1];
sx q[1];
rz(0.71539496) q[1];
x q[2];
rz(-0.69200055) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(-0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58671826) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2471152) q[0];
sx q[0];
rz(-3.101427) q[0];
sx q[0];
rz(-2.1263057) q[0];
rz(-2.9292551) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-0.36605787) q[1];
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
rz(-pi) q[1];
rz(1.4255849) q[2];
sx q[2];
rz(-1.8038857) q[2];
sx q[2];
rz(-0.22156246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8180994) q[1];
sx q[1];
rz(-2.202569) q[1];
sx q[1];
rz(1.9629823) q[1];
rz(-pi) q[2];
rz(2.8439431) q[3];
sx q[3];
rz(-2.1538056) q[3];
sx q[3];
rz(1.7109554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-0.69501957) q[2];
rz(0.85092893) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.2581576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224391) q[0];
sx q[0];
rz(-0.16534403) q[0];
sx q[0];
rz(-0.31546053) q[0];
rz(0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-0.42627898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75794719) q[0];
sx q[0];
rz(-1.8611127) q[0];
sx q[0];
rz(-2.2955672) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2832203) q[2];
sx q[2];
rz(-1.7933262) q[2];
sx q[2];
rz(-0.50642799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0002928) q[1];
sx q[1];
rz(-2.0054549) q[1];
sx q[1];
rz(-3.0463957) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8180646) q[3];
sx q[3];
rz(-2.0364722) q[3];
sx q[3];
rz(0.22899187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2686501) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(0.028701393) q[2];
rz(0.21102333) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(-0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(-2.9826214) q[0];
rz(2.4863415) q[1];
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
rz(2.4814674) q[0];
sx q[0];
rz(-0.9430389) q[0];
sx q[0];
rz(0.56076903) q[0];
x q[1];
rz(-1.9150429) q[2];
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
rz(-2.3895479) q[1];
sx q[1];
rz(-1.5524781) q[1];
sx q[1];
rz(2.1356029) q[1];
x q[2];
rz(-0.69055478) q[3];
sx q[3];
rz(-2.6752546) q[3];
sx q[3];
rz(-0.40377221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0752461) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(2.9998903) q[2];
rz(0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(-2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-0.53618479) q[0];
sx q[0];
rz(2.2313927) q[0];
rz(-2.3987112) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(1.9076617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155411) q[0];
sx q[0];
rz(-1.7103467) q[0];
sx q[0];
rz(1.6427137) q[0];
rz(-pi) q[1];
rz(-2.1338314) q[2];
sx q[2];
rz(-2.6673311) q[2];
sx q[2];
rz(-0.33111085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0764816) q[1];
sx q[1];
rz(-1.9533159) q[1];
sx q[1];
rz(0.72648813) q[1];
rz(-pi) q[2];
rz(-1.1804092) q[3];
sx q[3];
rz(-2.1258866) q[3];
sx q[3];
rz(-3.0267455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(0.91599715) q[2];
rz(-2.8902174) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6420355) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-0.04059759) q[0];
rz(-1.9675072) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(1.0601128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670676) q[0];
sx q[0];
rz(-0.39345783) q[0];
sx q[0];
rz(-2.8889546) q[0];
rz(1.4547075) q[2];
sx q[2];
rz(-1.5903859) q[2];
sx q[2];
rz(-1.3649114) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5225142) q[1];
sx q[1];
rz(-0.97303094) q[1];
sx q[1];
rz(1.1903711) q[1];
x q[2];
rz(-0.65877093) q[3];
sx q[3];
rz(-2.3537354) q[3];
sx q[3];
rz(0.066490563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4982831) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(-0.79116428) q[2];
rz(1.6061973) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(-0.78996381) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(2.9839363) q[0];
rz(0.12570307) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(-0.30473614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5253741) q[0];
sx q[0];
rz(-0.33712876) q[0];
sx q[0];
rz(1.824877) q[0];
x q[1];
rz(-0.53908277) q[2];
sx q[2];
rz(-0.7149019) q[2];
sx q[2];
rz(1.4366921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8915776) q[1];
sx q[1];
rz(-0.13430891) q[1];
sx q[1];
rz(3.1242643) q[1];
x q[2];
rz(2.7823506) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(-0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2711082) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.452272) q[2];
rz(-0.74635402) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8993503) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(-0.53264701) q[0];
rz(0.38231725) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(0.16758448) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2004799) q[0];
sx q[0];
rz(-0.22819209) q[0];
sx q[0];
rz(-2.1424293) q[0];
rz(2.8791076) q[2];
sx q[2];
rz(-2.1023316) q[2];
sx q[2];
rz(-2.558311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1780562) q[1];
sx q[1];
rz(-1.2884226) q[1];
sx q[1];
rz(2.8368188) q[1];
rz(-2.7277023) q[3];
sx q[3];
rz(-0.87855708) q[3];
sx q[3];
rz(-2.5138951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(-2.6746542) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(-1.9095437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(-2.3570428) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(-2.4541773) q[2];
sx q[2];
rz(-2.1599507) q[2];
sx q[2];
rz(0.31518515) q[2];
rz(-0.26950027) q[3];
sx q[3];
rz(-1.917358) q[3];
sx q[3];
rz(-2.6488398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
