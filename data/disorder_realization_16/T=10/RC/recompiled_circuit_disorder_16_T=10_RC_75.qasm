OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(0.33831236) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51988039) q[0];
sx q[0];
rz(-1.3268688) q[0];
sx q[0];
rz(1.8983311) q[0];
x q[1];
rz(2.2489684) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(-0.48722789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6602064) q[1];
sx q[1];
rz(-1.2848789) q[1];
sx q[1];
rz(1.5155161) q[1];
x q[2];
rz(2.0831574) q[3];
sx q[3];
rz(-1.5844987) q[3];
sx q[3];
rz(2.0126359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(-0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.8992791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466465) q[0];
sx q[0];
rz(-1.6992237) q[0];
sx q[0];
rz(2.8917679) q[0];
rz(-pi) q[1];
rz(-1.785948) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(1.3607963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17774432) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(-0.53667712) q[1];
rz(1.0817238) q[3];
sx q[3];
rz(-0.92182577) q[3];
sx q[3];
rz(-1.058941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(-0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-0.72845355) q[0];
rz(-1.6473673) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-2.1247991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8249614) q[0];
sx q[0];
rz(-1.5350071) q[0];
sx q[0];
rz(0.081598452) q[0];
rz(2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(1.0334894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9079202) q[1];
sx q[1];
rz(-2.3111812) q[1];
sx q[1];
rz(-2.9552712) q[1];
rz(3.0924762) q[3];
sx q[3];
rz(-1.2894863) q[3];
sx q[3];
rz(-2.1360872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(1.0495079) q[2];
rz(0.63878757) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(0.24681117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.4220337) q[0];
sx q[0];
rz(0.87417283) q[0];
x q[1];
rz(-2.6631782) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(-1.071655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2553195) q[1];
sx q[1];
rz(-1.6233994) q[1];
sx q[1];
rz(-2.7707151) q[1];
rz(-pi) q[2];
rz(-1.1770583) q[3];
sx q[3];
rz(-2.4871832) q[3];
sx q[3];
rz(0.27458336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(-1.8141618) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(-2.8932103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800956) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(-1.5310775) q[0];
rz(-pi) q[1];
rz(2.4400473) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(-1.0813576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3779113) q[1];
sx q[1];
rz(-1.7421725) q[1];
sx q[1];
rz(2.548449) q[1];
rz(0.6487209) q[3];
sx q[3];
rz(-1.3789163) q[3];
sx q[3];
rz(0.81024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(0.79745897) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(2.6388772) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1335063) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(0.1246917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42449441) q[0];
sx q[0];
rz(-1.4846804) q[0];
sx q[0];
rz(-2.9647102) q[0];
x q[1];
rz(2.5325534) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(-0.86953029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9166959) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(2.4042261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3521306) q[3];
sx q[3];
rz(-1.8018186) q[3];
sx q[3];
rz(-0.71393379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(2.52264) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-1.0707062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8585513) q[0];
sx q[0];
rz(-2.0548477) q[0];
sx q[0];
rz(-1.5701243) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7942065) q[2];
sx q[2];
rz(-1.6958478) q[2];
sx q[2];
rz(-1.8206247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5412022) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(-2.3669835) q[1];
rz(-pi) q[2];
rz(-1.5504863) q[3];
sx q[3];
rz(-1.2386525) q[3];
sx q[3];
rz(-0.16196812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038717) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(1.7659811) q[0];
x q[1];
rz(1.4119554) q[2];
sx q[2];
rz(-1.7454595) q[2];
sx q[2];
rz(-0.58089248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2944813) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(0.28114762) q[1];
rz(-pi) q[2];
rz(2.3341228) q[3];
sx q[3];
rz(-1.315457) q[3];
sx q[3];
rz(2.1047999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5802713) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(-2.8619213) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79265362) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(3.1025725) q[0];
rz(-pi) q[1];
rz(2.3896396) q[2];
sx q[2];
rz(-0.30138902) q[2];
sx q[2];
rz(-1.6981268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54535941) q[1];
sx q[1];
rz(-1.4701478) q[1];
sx q[1];
rz(-0.049493162) q[1];
rz(-pi) q[2];
rz(-0.46189724) q[3];
sx q[3];
rz(-1.4326722) q[3];
sx q[3];
rz(1.671333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(0.80741185) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40598665) q[0];
sx q[0];
rz(-1.1445023) q[0];
sx q[0];
rz(2.6398525) q[0];
x q[1];
rz(-0.21364613) q[2];
sx q[2];
rz(-1.7591957) q[2];
sx q[2];
rz(-1.633916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.909543) q[1];
sx q[1];
rz(-0.068040158) q[1];
sx q[1];
rz(1.0104936) q[1];
rz(-pi) q[2];
rz(-2.093408) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(-0.89695938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-0.87396809) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-0.59861029) q[2];
sx q[2];
rz(-0.90000464) q[2];
sx q[2];
rz(2.4760751) q[2];
rz(-2.279083) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
