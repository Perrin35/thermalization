OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(3.951374) q[0];
sx q[0];
rz(9.9561719) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(1.9385424) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8578313) q[0];
sx q[0];
rz(-0.28781048) q[0];
sx q[0];
rz(-2.3057111) q[0];
rz(1.3538023) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(-0.72460246) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6984451) q[1];
sx q[1];
rz(-1.6277715) q[1];
sx q[1];
rz(0.2351825) q[1];
x q[2];
rz(0.29294149) q[3];
sx q[3];
rz(-1.5353234) q[3];
sx q[3];
rz(0.5637416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.8751289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931912) q[0];
sx q[0];
rz(-1.9494434) q[0];
sx q[0];
rz(0.65403954) q[0];
rz(-pi) q[1];
rz(0.23898791) q[2];
sx q[2];
rz(-1.1766489) q[2];
sx q[2];
rz(0.1220526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5310865) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(-0.83341877) q[1];
x q[2];
rz(-2.0224138) q[3];
sx q[3];
rz(-2.0997117) q[3];
sx q[3];
rz(-0.12714735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1566029) q[2];
sx q[2];
rz(-1.9282324) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(1.7595694) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377797) q[0];
sx q[0];
rz(-0.63596361) q[0];
sx q[0];
rz(3.0657835) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1157672) q[2];
sx q[2];
rz(-2.2098944) q[2];
sx q[2];
rz(2.7884723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69933575) q[1];
sx q[1];
rz(-1.2067716) q[1];
sx q[1];
rz(2.060021) q[1];
rz(-pi) q[2];
rz(-2.2494499) q[3];
sx q[3];
rz(-0.8222848) q[3];
sx q[3];
rz(0.014978623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(1.5807318) q[2];
rz(2.1598699) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(-2.4981807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075539) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(-0.38726989) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(0.29104582) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211348) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(0.27340425) q[0];
rz(2.307345) q[2];
sx q[2];
rz(-0.4220037) q[2];
sx q[2];
rz(0.01854245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2070771) q[1];
sx q[1];
rz(-2.0743437) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(-pi) q[2];
rz(-2.980551) q[3];
sx q[3];
rz(-2.1412009) q[3];
sx q[3];
rz(-0.17810861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3691833) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(2.496526) q[0];
rz(-0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(2.5114139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5298313) q[0];
sx q[0];
rz(-2.1453619) q[0];
sx q[0];
rz(1.7163506) q[0];
rz(-pi) q[1];
rz(0.81996702) q[2];
sx q[2];
rz(-1.8659288) q[2];
sx q[2];
rz(-2.453014) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.899257) q[1];
sx q[1];
rz(-1.3471654) q[1];
sx q[1];
rz(-1.0244589) q[1];
rz(-2.583902) q[3];
sx q[3];
rz(-1.7002749) q[3];
sx q[3];
rz(1.492471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3877635) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(-1.2498614) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(2.2201339) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0548148) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(0.81370083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610634) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(-1.3570696) q[0];
x q[1];
rz(2.1505693) q[2];
sx q[2];
rz(-2.2057141) q[2];
sx q[2];
rz(2.6099176) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5972283) q[1];
sx q[1];
rz(-1.9042943) q[1];
sx q[1];
rz(-0.79798214) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86088647) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(-0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(2.4783321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2384757) q[0];
sx q[0];
rz(-0.65070063) q[0];
sx q[0];
rz(1.6102717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3528321) q[2];
sx q[2];
rz(-0.94467615) q[2];
sx q[2];
rz(1.2499714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85147334) q[1];
sx q[1];
rz(-2.0514601) q[1];
sx q[1];
rz(-2.9139247) q[1];
rz(-2.0459941) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(-1.0218395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(-0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-0.10738871) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-2.6838141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4360355) q[0];
sx q[0];
rz(-1.8920664) q[0];
sx q[0];
rz(0.9106439) q[0];
rz(-0.56577487) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(2.4772252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0486054) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(2.155817) q[1];
x q[2];
rz(-2.5910208) q[3];
sx q[3];
rz(-0.62567657) q[3];
sx q[3];
rz(-1.613137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.15437056) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-2.1881058) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-2.5568533) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027771587) q[0];
sx q[0];
rz(-2.6391533) q[0];
sx q[0];
rz(-2.9071945) q[0];
rz(-2.5768313) q[2];
sx q[2];
rz(-0.62099651) q[2];
sx q[2];
rz(-2.9850105) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.016420267) q[1];
sx q[1];
rz(-0.44508815) q[1];
sx q[1];
rz(-2.9801286) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14296602) q[3];
sx q[3];
rz(-1.8355832) q[3];
sx q[3];
rz(-0.20204443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-2.4861091) q[0];
rz(-1.1570702) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(1.5225333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64064202) q[0];
sx q[0];
rz(-0.92693936) q[0];
sx q[0];
rz(1.9615016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9117781) q[2];
sx q[2];
rz(-0.91170646) q[2];
sx q[2];
rz(1.1863208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5628221) q[1];
sx q[1];
rz(-2.4914503) q[1];
sx q[1];
rz(-0.27362089) q[1];
rz(-pi) q[2];
rz(2.7195815) q[3];
sx q[3];
rz(-1.2000298) q[3];
sx q[3];
rz(-1.5164204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15554252) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(2.299451) q[2];
rz(-1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-0.026528927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(2.4189667) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(-1.4459544) q[2];
sx q[2];
rz(-2.9308133) q[2];
sx q[2];
rz(1.4598893) q[2];
rz(0.079506569) q[3];
sx q[3];
rz(-0.88149298) q[3];
sx q[3];
rz(0.029824921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
