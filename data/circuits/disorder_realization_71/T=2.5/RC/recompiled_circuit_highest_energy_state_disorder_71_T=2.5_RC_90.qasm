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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(2.4655226) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(-0.91125429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9118496) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(0.020072083) q[0];
x q[1];
rz(0.010013117) q[2];
sx q[2];
rz(-0.62558936) q[2];
sx q[2];
rz(2.6715476) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.377502) q[1];
sx q[1];
rz(-1.5385748) q[1];
sx q[1];
rz(2.568214) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0082672) q[3];
sx q[3];
rz(-0.76078868) q[3];
sx q[3];
rz(-2.4530792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(2.3294219) q[2];
rz(-0.075856097) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(2.5465452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512063) q[0];
sx q[0];
rz(-0.47845978) q[0];
sx q[0];
rz(1.867021) q[0];
rz(1.5050585) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(1.1638181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75117373) q[0];
sx q[0];
rz(-1.112022) q[0];
sx q[0];
rz(-1.0659046) q[0];
rz(-2.4435448) q[2];
sx q[2];
rz(-2.5155332) q[2];
sx q[2];
rz(1.9174089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5806942) q[1];
sx q[1];
rz(-2.1836062) q[1];
sx q[1];
rz(1.3152907) q[1];
rz(-pi) q[2];
rz(-1.9979565) q[3];
sx q[3];
rz(-1.2795951) q[3];
sx q[3];
rz(1.8607651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.290648) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5811496) q[2];
rz(-2.9110939) q[3];
sx q[3];
rz(-2.0932308) q[3];
sx q[3];
rz(-2.7010664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59349638) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(3.0189959) q[0];
rz(-2.3444046) q[1];
sx q[1];
rz(-2.1295348) q[1];
sx q[1];
rz(-3.0534993) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0038550475) q[0];
sx q[0];
rz(-0.81043506) q[0];
sx q[0];
rz(1.9268981) q[0];
rz(-pi) q[1];
rz(0.083115301) q[2];
sx q[2];
rz(-1.7975472) q[2];
sx q[2];
rz(-0.47435681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5777307) q[1];
sx q[1];
rz(-0.20521611) q[1];
sx q[1];
rz(-0.83744253) q[1];
rz(-pi) q[2];
rz(-2.5949305) q[3];
sx q[3];
rz(-2.529749) q[3];
sx q[3];
rz(-0.33239588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2587571) q[2];
sx q[2];
rz(-1.5587403) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(2.6435408) q[3];
sx q[3];
rz(-0.96314722) q[3];
sx q[3];
rz(-2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8525456) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(2.9943941) q[0];
rz(2.4920801) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(1.46896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8627074) q[0];
sx q[0];
rz(-1.5636347) q[0];
sx q[0];
rz(-1.7938406) q[0];
rz(-pi) q[1];
rz(1.1873522) q[2];
sx q[2];
rz(-0.45971515) q[2];
sx q[2];
rz(-1.6123259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4893579) q[1];
sx q[1];
rz(-1.2235886) q[1];
sx q[1];
rz(-0.90851291) q[1];
rz(-0.62375237) q[3];
sx q[3];
rz(-0.68188462) q[3];
sx q[3];
rz(-2.9870913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2088251) q[2];
sx q[2];
rz(-2.6471477) q[2];
sx q[2];
rz(0.24173582) q[2];
rz(-0.73511165) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(2.2562476) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(-0.56040323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86515988) q[0];
sx q[0];
rz(-1.393228) q[0];
sx q[0];
rz(-1.0516758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8444195) q[2];
sx q[2];
rz(-1.3577596) q[2];
sx q[2];
rz(1.7293255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.4679298) q[1];
sx q[1];
rz(-2.0286252) q[1];
sx q[1];
rz(-0.28854847) q[1];
x q[2];
rz(0.82920977) q[3];
sx q[3];
rz(-2.5253339) q[3];
sx q[3];
rz(2.5065638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6178599) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.3699646) q[3];
sx q[3];
rz(3.052616) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646249) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(-2.4846039) q[0];
rz(-0.14343801) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(0.73854804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8764832) q[0];
sx q[0];
rz(-0.68317693) q[0];
sx q[0];
rz(-2.9580412) q[0];
rz(2.9159268) q[2];
sx q[2];
rz(-1.0899001) q[2];
sx q[2];
rz(1.7250329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1164897) q[1];
sx q[1];
rz(-1.8855576) q[1];
sx q[1];
rz(0.88309755) q[1];
rz(-pi) q[2];
rz(0.62885999) q[3];
sx q[3];
rz(-2.326528) q[3];
sx q[3];
rz(0.12303837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8831545) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(0.0030041791) q[2];
rz(-0.7891807) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(1.167231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060519144) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(-0.14807598) q[0];
rz(-0.5332467) q[1];
sx q[1];
rz(-0.24594578) q[1];
sx q[1];
rz(-1.6291133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0740967) q[0];
sx q[0];
rz(-1.8379015) q[0];
sx q[0];
rz(2.6617034) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.327269) q[2];
sx q[2];
rz(-1.0934798) q[2];
sx q[2];
rz(-2.5392591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.032669205) q[1];
sx q[1];
rz(-1.2315688) q[1];
sx q[1];
rz(-3.0134588) q[1];
rz(-0.21085994) q[3];
sx q[3];
rz(-0.85975826) q[3];
sx q[3];
rz(-0.27121997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6988354) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(3.0485145) q[2];
rz(-0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(-1.9183581) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020697866) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(0.66463941) q[0];
rz(1.3497893) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(-0.68914366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1375547) q[0];
sx q[0];
rz(-1.327164) q[0];
sx q[0];
rz(1.6904657) q[0];
rz(2.1470239) q[2];
sx q[2];
rz(-2.3625018) q[2];
sx q[2];
rz(1.7528723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41708159) q[1];
sx q[1];
rz(-1.8838091) q[1];
sx q[1];
rz(-1.4177807) q[1];
x q[2];
rz(-0.34948055) q[3];
sx q[3];
rz(-2.9664024) q[3];
sx q[3];
rz(-0.82621208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.050921116) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(-1.3760759) q[2];
rz(-1.225166) q[3];
sx q[3];
rz(-2.4296032) q[3];
sx q[3];
rz(3.0998949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10251481) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(0.59120375) q[0];
rz(-2.8996331) q[1];
sx q[1];
rz(-2.0457025) q[1];
sx q[1];
rz(-2.718149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5862783) q[0];
sx q[0];
rz(-1.269377) q[0];
sx q[0];
rz(2.8239646) q[0];
rz(0.19337635) q[2];
sx q[2];
rz(-0.49623734) q[2];
sx q[2];
rz(-1.9691182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6953227) q[1];
sx q[1];
rz(-2.5188418) q[1];
sx q[1];
rz(2.0957309) q[1];
rz(-0.61708151) q[3];
sx q[3];
rz(-1.2344196) q[3];
sx q[3];
rz(-1.468889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5847136) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(1.1663743) q[2];
rz(1.1276468) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(-2.6166272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135076) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(-2.4684913) q[0];
rz(0.85064864) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(-2.8180715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44064012) q[0];
sx q[0];
rz(-1.8201168) q[0];
sx q[0];
rz(-1.9560157) q[0];
x q[1];
rz(-2.0704248) q[2];
sx q[2];
rz(-0.53348225) q[2];
sx q[2];
rz(0.96772742) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0464152) q[1];
sx q[1];
rz(-0.84734619) q[1];
sx q[1];
rz(-2.0207538) q[1];
x q[2];
rz(1.5629285) q[3];
sx q[3];
rz(-2.1546225) q[3];
sx q[3];
rz(1.5567832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(-0.072048135) q[2];
rz(2.1273023) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2608248) q[0];
sx q[0];
rz(-1.7171971) q[0];
sx q[0];
rz(1.1203753) q[0];
rz(2.8063759) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(-2.0497132) q[2];
sx q[2];
rz(-1.6960232) q[2];
sx q[2];
rz(-2.1868119) q[2];
rz(2.9670197) q[3];
sx q[3];
rz(-2.9963507) q[3];
sx q[3];
rz(0.33490845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
