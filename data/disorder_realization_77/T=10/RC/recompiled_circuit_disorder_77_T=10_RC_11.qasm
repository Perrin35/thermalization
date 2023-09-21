OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-2.9021184) q[0];
x q[1];
rz(0.92476966) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(-0.31121635) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9819298) q[1];
sx q[1];
rz(-2.5291981) q[1];
sx q[1];
rz(-2.0484522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16529103) q[3];
sx q[3];
rz(-2.878302) q[3];
sx q[3];
rz(0.96112448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1412927) q[0];
sx q[0];
rz(-0.89501689) q[0];
sx q[0];
rz(2.7594901) q[0];
rz(-0.58354124) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5838503) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18560219) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(-0.76693265) q[1];
rz(0.50206708) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(-0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-0.056578606) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.409257) q[0];
sx q[0];
rz(-2.0154672) q[0];
sx q[0];
rz(-0.94629855) q[0];
rz(-pi) q[1];
rz(0.25643202) q[2];
sx q[2];
rz(-1.7000546) q[2];
sx q[2];
rz(-1.749922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1287071) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(2.835564) q[1];
rz(-1.1832841) q[3];
sx q[3];
rz(-0.96252493) q[3];
sx q[3];
rz(-1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(-1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0915506) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(2.0043623) q[0];
x q[1];
rz(1.6323339) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(0.57519826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5220118) q[1];
sx q[1];
rz(-1.3995692) q[1];
sx q[1];
rz(-0.79230688) q[1];
rz(-2.2673244) q[3];
sx q[3];
rz(-2.4618751) q[3];
sx q[3];
rz(3.0391039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-0.11894225) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34236318) q[0];
sx q[0];
rz(-0.21393299) q[0];
sx q[0];
rz(-1.7844723) q[0];
x q[1];
rz(0.40446754) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(-2.5596465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7625092) q[1];
sx q[1];
rz(-3.0804539) q[1];
sx q[1];
rz(1.5361384) q[1];
x q[2];
rz(-0.43290187) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-0.054919682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024825) q[0];
sx q[0];
rz(-1.6533028) q[0];
sx q[0];
rz(-1.301618) q[0];
rz(-pi) q[1];
rz(-1.8940077) q[2];
sx q[2];
rz(-1.5885457) q[2];
sx q[2];
rz(-0.36349597) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1153032) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(-2.5251212) q[1];
rz(-2.218607) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(-1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(0.91032666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658265) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.5058917) q[0];
x q[1];
rz(0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(0.90781462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(-1.7757225) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4207553) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(0.48203326) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3807316) q[2];
sx q[2];
rz(-1.0957452) q[2];
sx q[2];
rz(0.27981112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(-0.47972958) q[1];
x q[2];
rz(-1.9037876) q[3];
sx q[3];
rz(-2.3559642) q[3];
sx q[3];
rz(1.1803407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(0.78197455) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-0.75884563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.471506) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(1.612081) q[0];
rz(-pi) q[1];
rz(-1.3779638) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(-2.6892975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5938877) q[1];
sx q[1];
rz(-2.2624359) q[1];
sx q[1];
rz(1.3681075) q[1];
rz(-pi) q[2];
rz(3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-0.60992253) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(-2.7479991) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2953193) q[2];
sx q[2];
rz(-1.2071949) q[2];
sx q[2];
rz(-1.3678577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6590609) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-0.34300967) q[1];
rz(-pi) q[2];
rz(1.0697332) q[3];
sx q[3];
rz(-2.2483629) q[3];
sx q[3];
rz(-0.60615218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80355766) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(2.4907885) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.4177633) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
rz(-2.0416904) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
