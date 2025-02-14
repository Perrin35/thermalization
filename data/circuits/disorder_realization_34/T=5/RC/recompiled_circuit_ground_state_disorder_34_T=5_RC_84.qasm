OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(-2.5791383) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0732361) q[0];
sx q[0];
rz(-2.773065) q[0];
sx q[0];
rz(1.4224398) q[0];
x q[1];
rz(-1.063827) q[2];
sx q[2];
rz(-1.7071144) q[2];
sx q[2];
rz(-0.16198128) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48784287) q[1];
sx q[1];
rz(-1.1636793) q[1];
sx q[1];
rz(-0.458325) q[1];
rz(1.1645557) q[3];
sx q[3];
rz(-1.2066168) q[3];
sx q[3];
rz(-1.4573569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(-0.79929024) q[2];
rz(-0.47131395) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(-1.348173) q[0];
rz(1.7680291) q[1];
sx q[1];
rz(-1.9919688) q[1];
sx q[1];
rz(1.9893533) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2179759) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(-2.7569735) q[0];
x q[1];
rz(0.94448467) q[2];
sx q[2];
rz(-0.99530333) q[2];
sx q[2];
rz(1.6863509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5785529) q[1];
sx q[1];
rz(-2.4336947) q[1];
sx q[1];
rz(2.8824174) q[1];
x q[2];
rz(-2.6226085) q[3];
sx q[3];
rz(-2.0618084) q[3];
sx q[3];
rz(1.6512914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(0.10654199) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(-2.54134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057864144) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-2.6237543) q[0];
rz(-0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(-2.6944366) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8661506) q[0];
sx q[0];
rz(-0.26922154) q[0];
sx q[0];
rz(0.24193995) q[0];
x q[1];
rz(-2.1250399) q[2];
sx q[2];
rz(-0.92281658) q[2];
sx q[2];
rz(-2.3584443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4616926) q[1];
sx q[1];
rz(-0.60631207) q[1];
sx q[1];
rz(1.3316505) q[1];
rz(-pi) q[2];
rz(-2.6833862) q[3];
sx q[3];
rz(-2.5669006) q[3];
sx q[3];
rz(1.1632869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7670224) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(2.6089597) q[2];
rz(0.24752188) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(-2.1029162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5175051) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(-3.0349773) q[0];
rz(-2.7952349) q[1];
sx q[1];
rz(-0.84500161) q[1];
sx q[1];
rz(1.1603629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713523) q[0];
sx q[0];
rz(-1.3246857) q[0];
sx q[0];
rz(1.7497803) q[0];
x q[1];
rz(2.423942) q[2];
sx q[2];
rz(-2.1776878) q[2];
sx q[2];
rz(0.90536149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6667156) q[1];
sx q[1];
rz(-1.8711539) q[1];
sx q[1];
rz(1.9431861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0456234) q[3];
sx q[3];
rz(-1.1743288) q[3];
sx q[3];
rz(-0.1016271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13566636) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(-1.7990254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0214486) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(-2.7914877) q[0];
rz(-0.95589751) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.9116481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1010001) q[0];
sx q[0];
rz(-0.44389566) q[0];
sx q[0];
rz(2.1568265) q[0];
x q[1];
rz(-1.4054662) q[2];
sx q[2];
rz(-1.5348736) q[2];
sx q[2];
rz(-1.3816116) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6767439) q[1];
sx q[1];
rz(-2.7584761) q[1];
sx q[1];
rz(-2.0967469) q[1];
rz(-pi) q[2];
rz(0.98389174) q[3];
sx q[3];
rz(-1.7376889) q[3];
sx q[3];
rz(1.9092164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(-2.7367075) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(-2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006007) q[0];
sx q[0];
rz(-1.0772935) q[0];
sx q[0];
rz(-0.58240044) q[0];
rz(-0.68663418) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(-1.2581717) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95520407) q[0];
sx q[0];
rz(-1.135728) q[0];
sx q[0];
rz(0.13829903) q[0];
rz(2.1297087) q[2];
sx q[2];
rz(-2.5311573) q[2];
sx q[2];
rz(1.287078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1327136) q[1];
sx q[1];
rz(-3.0522484) q[1];
sx q[1];
rz(2.4650399) q[1];
x q[2];
rz(0.21922501) q[3];
sx q[3];
rz(-0.98837438) q[3];
sx q[3];
rz(-1.5150013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(1.0180391) q[2];
rz(2.8954519) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(-0.14353453) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-2.9339583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.62791) q[0];
sx q[0];
rz(-1.915689) q[0];
sx q[0];
rz(2.3570127) q[0];
x q[1];
rz(0.34822627) q[2];
sx q[2];
rz(-2.1219606) q[2];
sx q[2];
rz(1.7334235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.138306) q[1];
sx q[1];
rz(-0.7423183) q[1];
sx q[1];
rz(-1.9600541) q[1];
rz(-pi) q[2];
rz(0.36334857) q[3];
sx q[3];
rz(-1.4814113) q[3];
sx q[3];
rz(2.9575728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2988854) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(-0.60619727) q[2];
rz(-0.68743622) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48924482) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(1.0580753) q[0];
rz(0.10969133) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(-1.6995957) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82175151) q[0];
sx q[0];
rz(-1.1761001) q[0];
sx q[0];
rz(-1.9072397) q[0];
x q[1];
rz(0.47693129) q[2];
sx q[2];
rz(-0.55185807) q[2];
sx q[2];
rz(-2.1546006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86275872) q[1];
sx q[1];
rz(-0.93288427) q[1];
sx q[1];
rz(0.61418038) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2992925) q[3];
sx q[3];
rz(-1.3437004) q[3];
sx q[3];
rz(-0.61448594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.480392) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(1.2083758) q[2];
rz(-2.4783573) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709568) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(0.092967689) q[0];
rz(-1.2804821) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(0.13519898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2988578) q[0];
sx q[0];
rz(-1.6162685) q[0];
sx q[0];
rz(3.088515) q[0];
rz(-0.5298631) q[2];
sx q[2];
rz(-2.8447897) q[2];
sx q[2];
rz(2.5713845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1982806) q[1];
sx q[1];
rz(-1.824728) q[1];
sx q[1];
rz(-1.5928245) q[1];
x q[2];
rz(0.074196176) q[3];
sx q[3];
rz(-1.7230986) q[3];
sx q[3];
rz(-2.0366813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.9229527) q[2];
sx q[2];
rz(-2.3700628) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58525697) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(-1.1248032) q[0];
rz(-1.8428165) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(-2.3769456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9287531) q[0];
sx q[0];
rz(-3.0055586) q[0];
sx q[0];
rz(-0.68599756) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49955274) q[2];
sx q[2];
rz(-2.1883114) q[2];
sx q[2];
rz(-0.90047405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5717126) q[1];
sx q[1];
rz(-0.24001828) q[1];
sx q[1];
rz(-1.0102788) q[1];
rz(-pi) q[2];
rz(-1.1922791) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(-1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3760959) q[2];
sx q[2];
rz(-1.2906047) q[2];
sx q[2];
rz(2.9343904) q[2];
rz(2.1616705) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021066396) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(-1.4516713) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(-0.23811447) q[3];
sx q[3];
rz(-1.7287711) q[3];
sx q[3];
rz(0.14915376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
