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
rz(-0.67607003) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(-0.91125429) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.229743) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(0.020072083) q[0];
rz(-pi) q[1];
rz(0.62556556) q[2];
sx q[2];
rz(-1.5766597) q[2];
sx q[2];
rz(2.0489583) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76409066) q[1];
sx q[1];
rz(-1.6030178) q[1];
sx q[1];
rz(-2.568214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6966713) q[3];
sx q[3];
rz(-2.3231703) q[3];
sx q[3];
rz(-0.5054121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91997826) q[2];
sx q[2];
rz(-0.45404926) q[2];
sx q[2];
rz(-2.3294219) q[2];
rz(0.075856097) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(-2.5465452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6512063) q[0];
sx q[0];
rz(-0.47845978) q[0];
sx q[0];
rz(1.2745717) q[0];
rz(-1.5050585) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(-1.1638181) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620016) q[0];
sx q[0];
rz(-1.1222071) q[0];
sx q[0];
rz(2.6278087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6357119) q[2];
sx q[2];
rz(-1.1846647) q[2];
sx q[2];
rz(-0.94368151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5608985) q[1];
sx q[1];
rz(-2.1836062) q[1];
sx q[1];
rz(-1.3152907) q[1];
x q[2];
rz(-2.8234661) q[3];
sx q[3];
rz(-1.162718) q[3];
sx q[3];
rz(-2.7216743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.290648) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5811496) q[2];
rz(0.23049878) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(2.7010664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59349638) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(0.12259677) q[0];
rz(-2.3444046) q[1];
sx q[1];
rz(-2.1295348) q[1];
sx q[1];
rz(0.0880934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255955) q[0];
sx q[0];
rz(-1.8261709) q[0];
sx q[0];
rz(2.3488464) q[0];
x q[1];
rz(0.083115301) q[2];
sx q[2];
rz(-1.3440455) q[2];
sx q[2];
rz(-2.6672358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3077724) q[1];
sx q[1];
rz(-1.4188179) q[1];
sx q[1];
rz(-0.13843468) q[1];
rz(2.6016079) q[3];
sx q[3];
rz(-1.8740046) q[3];
sx q[3];
rz(1.4410415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(-1.8176414) q[2];
rz(-2.6435408) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(-2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8525456) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(2.9943941) q[0];
rz(0.64951253) q[1];
sx q[1];
rz(-1.6308866) q[1];
sx q[1];
rz(1.46896) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32347476) q[0];
sx q[0];
rz(-2.9184353) q[0];
sx q[0];
rz(-1.5384307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.001226) q[2];
sx q[2];
rz(-1.4040315) q[2];
sx q[2];
rz(2.8362136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33002871) q[1];
sx q[1];
rz(-0.73544502) q[1];
sx q[1];
rz(2.1027019) q[1];
rz(-2.0135512) q[3];
sx q[3];
rz(-2.1078133) q[3];
sx q[3];
rz(-2.2396416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9327675) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(-2.8998568) q[2];
rz(2.406481) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46297896) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(3.0188634) q[0];
rz(-0.8853451) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(-0.56040323) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80622506) q[0];
sx q[0];
rz(-1.0606375) q[0];
sx q[0];
rz(-2.9377743) q[0];
x q[1];
rz(0.29717314) q[2];
sx q[2];
rz(-1.783833) q[2];
sx q[2];
rz(-1.7293255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1691815) q[1];
sx q[1];
rz(-1.3126846) q[1];
sx q[1];
rz(-1.0960137) q[1];
rz(-pi) q[2];
rz(1.089483) q[3];
sx q[3];
rz(-1.9718687) q[3];
sx q[3];
rz(0.29386753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5237328) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8646249) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(-0.65698874) q[0];
rz(-0.14343801) q[1];
sx q[1];
rz(-1.5851494) q[1];
sx q[1];
rz(2.4030446) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114233) q[0];
sx q[0];
rz(-2.2403754) q[0];
sx q[0];
rz(-1.4233146) q[0];
x q[1];
rz(0.22566585) q[2];
sx q[2];
rz(-1.0899001) q[2];
sx q[2];
rz(1.4165598) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.025102928) q[1];
sx q[1];
rz(-1.2560351) q[1];
sx q[1];
rz(-2.2584951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62885999) q[3];
sx q[3];
rz(-2.326528) q[3];
sx q[3];
rz(0.12303837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2584381) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(-3.1385885) q[2];
rz(2.352412) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(-1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0810735) q[0];
sx q[0];
rz(-2.7923212) q[0];
sx q[0];
rz(2.9935167) q[0];
rz(-0.5332467) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(1.6291133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142198) q[0];
sx q[0];
rz(-0.54410579) q[0];
sx q[0];
rz(2.6065488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.523411) q[2];
sx q[2];
rz(-2.2264109) q[2];
sx q[2];
rz(0.55936343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.032669205) q[1];
sx q[1];
rz(-1.2315688) q[1];
sx q[1];
rz(-3.0134588) q[1];
x q[2];
rz(-0.21085994) q[3];
sx q[3];
rz(-2.2818344) q[3];
sx q[3];
rz(-2.8703727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44275722) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(0.093078144) q[2];
rz(0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.9183581) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020697866) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(-0.66463941) q[0];
rz(1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(0.68914366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5377602) q[0];
sx q[0];
rz(-1.6869154) q[0];
sx q[0];
rz(0.24531792) q[0];
rz(0.99456878) q[2];
sx q[2];
rz(-2.3625018) q[2];
sx q[2];
rz(-1.7528723) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0943567) q[1];
sx q[1];
rz(-2.7942889) q[1];
sx q[1];
rz(-0.44012614) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1647968) q[3];
sx q[3];
rz(-1.6305123) q[3];
sx q[3];
rz(-2.7415558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0906715) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(1.3760759) q[2];
rz(-1.225166) q[3];
sx q[3];
rz(-2.4296032) q[3];
sx q[3];
rz(-0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390778) q[0];
sx q[0];
rz(-3.097105) q[0];
sx q[0];
rz(-0.59120375) q[0];
rz(-0.2419596) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(-2.718149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11277448) q[0];
sx q[0];
rz(-1.2679546) q[0];
sx q[0];
rz(1.2545198) q[0];
rz(-pi) q[1];
rz(-2.9482163) q[2];
sx q[2];
rz(-0.49623734) q[2];
sx q[2];
rz(1.1724745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.075942) q[1];
sx q[1];
rz(-1.0417176) q[1];
sx q[1];
rz(0.34543646) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9758281) q[3];
sx q[3];
rz(-0.99294239) q[3];
sx q[3];
rz(-0.33194968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55687904) q[2];
sx q[2];
rz(-0.60606474) q[2];
sx q[2];
rz(1.9752183) q[2];
rz(2.0139458) q[3];
sx q[3];
rz(-0.96846628) q[3];
sx q[3];
rz(-2.6166272) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128085) q[0];
sx q[0];
rz(-2.7078244) q[0];
sx q[0];
rz(-2.4684913) q[0];
rz(-2.290944) q[1];
sx q[1];
rz(-1.7885845) q[1];
sx q[1];
rz(2.8180715) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1111543) q[0];
sx q[0];
rz(-1.9435099) q[0];
sx q[0];
rz(-0.26813676) q[0];
x q[1];
rz(0.27574972) q[2];
sx q[2];
rz(-1.1080896) q[2];
sx q[2];
rz(-1.6088865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4648351) q[1];
sx q[1];
rz(-0.82993648) q[1];
sx q[1];
rz(2.6840032) q[1];
x q[2];
rz(1.5786641) q[3];
sx q[3];
rz(-0.98697013) q[3];
sx q[3];
rz(1.5567832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2181776) q[2];
sx q[2];
rz(-2.9604993) q[2];
sx q[2];
rz(0.072048135) q[2];
rz(-1.0142903) q[3];
sx q[3];
rz(-0.91599661) q[3];
sx q[3];
rz(0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(1.3041244) q[2];
sx q[2];
rz(-2.647807) q[2];
sx q[2];
rz(-0.37995445) q[2];
rz(1.5961965) q[3];
sx q[3];
rz(-1.7138154) q[3];
sx q[3];
rz(-2.9830767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
