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
rz(0.060184181) q[0];
sx q[0];
rz(-2.0826075) q[0];
sx q[0];
rz(-1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84278569) q[0];
sx q[0];
rz(-1.1514542) q[0];
sx q[0];
rz(1.4489277) q[0];
rz(-pi) q[1];
rz(-3.0344738) q[2];
sx q[2];
rz(-0.20640443) q[2];
sx q[2];
rz(2.602488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0204326) q[1];
sx q[1];
rz(-1.766664) q[1];
sx q[1];
rz(-1.905026) q[1];
x q[2];
rz(-0.20186353) q[3];
sx q[3];
rz(-2.2998357) q[3];
sx q[3];
rz(2.1191459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4787204) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(-1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(-0.61479968) q[0];
rz(-2.2531033) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(2.6420171) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0083778) q[0];
sx q[0];
rz(-1.7724123) q[0];
sx q[0];
rz(1.9033543) q[0];
rz(-2.8728054) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(-0.42551431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0608637) q[1];
sx q[1];
rz(-1.7434967) q[1];
sx q[1];
rz(-2.2279146) q[1];
rz(-pi) q[2];
rz(2.0604443) q[3];
sx q[3];
rz(-0.64125618) q[3];
sx q[3];
rz(1.412751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(0.020817967) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-3.0121682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11331081) q[0];
sx q[0];
rz(-0.89963642) q[0];
sx q[0];
rz(-2.9379803) q[0];
x q[1];
rz(2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(-2.3455623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30363917) q[1];
sx q[1];
rz(-0.15175125) q[1];
sx q[1];
rz(0.93317731) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.09472) q[3];
sx q[3];
rz(-2.2942846) q[3];
sx q[3];
rz(2.0499961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0176598) q[2];
sx q[2];
rz(-1.5256226) q[2];
sx q[2];
rz(-1.370149) q[2];
rz(2.395199) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(-0.56513894) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(2.3133004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6713326) q[0];
sx q[0];
rz(-1.4909706) q[0];
sx q[0];
rz(0.74175394) q[0];
rz(-1.0617375) q[2];
sx q[2];
rz(-2.5974496) q[2];
sx q[2];
rz(-0.5296112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9984549) q[1];
sx q[1];
rz(-1.3052485) q[1];
sx q[1];
rz(2.6396855) q[1];
rz(-1.9902112) q[3];
sx q[3];
rz(-0.55636251) q[3];
sx q[3];
rz(-1.1384979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(-1.431541) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(-1.8121207) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(0.00024814127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.3868185) q[0];
sx q[0];
rz(-1.8198245) q[0];
rz(-pi) q[1];
rz(-3.0543112) q[2];
sx q[2];
rz(-1.0382663) q[2];
sx q[2];
rz(-2.0244348) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8793638) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(-1.0443347) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4100634) q[3];
sx q[3];
rz(-1.6247182) q[3];
sx q[3];
rz(1.9526854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1382711) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936845) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(-0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(2.1991275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346325) q[0];
sx q[0];
rz(-0.31123268) q[0];
sx q[0];
rz(-0.65629058) q[0];
rz(-pi) q[1];
rz(-1.8293657) q[2];
sx q[2];
rz(-1.2814008) q[2];
sx q[2];
rz(-1.0602151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2561431) q[1];
sx q[1];
rz(-1.1032618) q[1];
sx q[1];
rz(2.8199151) q[1];
rz(0.26047892) q[3];
sx q[3];
rz(-2.0325869) q[3];
sx q[3];
rz(2.8094629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59948644) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.2055509) q[3];
sx q[3];
rz(2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(3.0410774) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(2.2844792) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53062886) q[0];
sx q[0];
rz(-0.93994323) q[0];
sx q[0];
rz(1.1683589) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4947462) q[2];
sx q[2];
rz(-1.5313206) q[2];
sx q[2];
rz(-1.0919065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3639314) q[1];
sx q[1];
rz(-0.45501935) q[1];
sx q[1];
rz(-0.48952405) q[1];
rz(1.1825652) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.7804467) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(-0.030990344) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.2773638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63420682) q[0];
sx q[0];
rz(-1.2370101) q[0];
sx q[0];
rz(2.7884363) q[0];
rz(-pi) q[1];
rz(-2.8988016) q[2];
sx q[2];
rz(-2.0567679) q[2];
sx q[2];
rz(2.9832341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.053169202) q[1];
sx q[1];
rz(-0.49234566) q[1];
sx q[1];
rz(1.4150934) q[1];
x q[2];
rz(-1.8559009) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(1.7948732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.10451) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(-2.6541397) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(1.4016271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6232672) q[0];
sx q[0];
rz(-1.6931769) q[0];
sx q[0];
rz(-2.8920679) q[0];
rz(2.6675111) q[2];
sx q[2];
rz(-2.3326121) q[2];
sx q[2];
rz(-0.37341973) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9822787) q[1];
sx q[1];
rz(-1.754292) q[1];
sx q[1];
rz(-0.029551701) q[1];
rz(-pi) q[2];
rz(2.1622439) q[3];
sx q[3];
rz(-0.92100793) q[3];
sx q[3];
rz(1.5486919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(0.6366716) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.3670032) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(-1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.6148753) q[1];
sx q[1];
rz(1.0640594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2702613) q[0];
sx q[0];
rz(-0.0055905213) q[0];
sx q[0];
rz(-2.5020775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80697717) q[2];
sx q[2];
rz(-0.39013559) q[2];
sx q[2];
rz(1.5240508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56434332) q[1];
sx q[1];
rz(-1.6112304) q[1];
sx q[1];
rz(-3.0811429) q[1];
x q[2];
rz(1.6056772) q[3];
sx q[3];
rz(-1.1115171) q[3];
sx q[3];
rz(1.9495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4260063) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(-3.0675724) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1114125) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(0.61182712) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(0.62803531) q[2];
sx q[2];
rz(-0.59878329) q[2];
sx q[2];
rz(-1.5046635) q[2];
rz(-1.0985804) q[3];
sx q[3];
rz(-2.2618812) q[3];
sx q[3];
rz(-1.19899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
