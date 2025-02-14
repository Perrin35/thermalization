OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3235029) q[0];
sx q[0];
rz(-2.7348195) q[0];
sx q[0];
rz(-0.50954252) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(6.4767467) q[1];
sx q[1];
rz(11.560796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60765841) q[0];
sx q[0];
rz(-0.21398057) q[0];
sx q[0];
rz(2.7817552) q[0];
x q[1];
rz(-0.92410134) q[2];
sx q[2];
rz(-1.8515203) q[2];
sx q[2];
rz(-1.678987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79709681) q[1];
sx q[1];
rz(-2.2386947) q[1];
sx q[1];
rz(2.383166) q[1];
rz(-pi) q[2];
rz(2.3629684) q[3];
sx q[3];
rz(-1.1168241) q[3];
sx q[3];
rz(-3.0954645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84844184) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(2.0330009) q[2];
rz(-2.0075924) q[3];
sx q[3];
rz(-2.0847376) q[3];
sx q[3];
rz(-3.0602684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14652458) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(2.8296237) q[0];
rz(-2.9106855) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(-1.1280967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5018092) q[0];
sx q[0];
rz(-2.6310496) q[0];
sx q[0];
rz(-1.9656409) q[0];
rz(2.6083469) q[2];
sx q[2];
rz(-0.74138481) q[2];
sx q[2];
rz(-2.7864151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.038874168) q[1];
sx q[1];
rz(-2.6800214) q[1];
sx q[1];
rz(-0.25074236) q[1];
rz(-2.6873782) q[3];
sx q[3];
rz(-0.81791211) q[3];
sx q[3];
rz(2.1061768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4565304) q[2];
sx q[2];
rz(-1.5038467) q[2];
sx q[2];
rz(-0.064229639) q[2];
rz(-1.8188933) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(-2.2548811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(1.892426) q[0];
rz(2.8104172) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(1.9452728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099469378) q[0];
sx q[0];
rz(-1.5705918) q[0];
sx q[0];
rz(-1.571079) q[0];
rz(-0.021042391) q[2];
sx q[2];
rz(-2.8940563) q[2];
sx q[2];
rz(2.4602082) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0374245) q[1];
sx q[1];
rz(-2.0817887) q[1];
sx q[1];
rz(-0.68242208) q[1];
rz(1.182797) q[3];
sx q[3];
rz(-0.85577589) q[3];
sx q[3];
rz(-1.9388461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(-0.91274846) q[2];
rz(-1.7031472) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(-1.8057757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790633) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(0.67034876) q[0];
rz(-0.66728512) q[1];
sx q[1];
rz(-1.0083464) q[1];
sx q[1];
rz(-0.60428062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9491228) q[0];
sx q[0];
rz(-1.5274824) q[0];
sx q[0];
rz(-2.4049824) q[0];
rz(-2.759614) q[2];
sx q[2];
rz(-1.7151217) q[2];
sx q[2];
rz(1.195418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5126219) q[1];
sx q[1];
rz(-1.781946) q[1];
sx q[1];
rz(1.8308012) q[1];
rz(2.3015062) q[3];
sx q[3];
rz(-1.4021676) q[3];
sx q[3];
rz(-0.11883277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1112572) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(2.7643909) q[2];
rz(0.12353573) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(1.7326573) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971624) q[0];
sx q[0];
rz(-0.57752174) q[0];
sx q[0];
rz(-2.8422728) q[0];
rz(-2.0206644) q[1];
sx q[1];
rz(-1.2052373) q[1];
sx q[1];
rz(-2.1790806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045557307) q[0];
sx q[0];
rz(-1.0518414) q[0];
sx q[0];
rz(-0.23764289) q[0];
rz(0.91470529) q[2];
sx q[2];
rz(-1.1924628) q[2];
sx q[2];
rz(1.8112184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1735517) q[1];
sx q[1];
rz(-2.3690201) q[1];
sx q[1];
rz(-0.79141683) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5425517) q[3];
sx q[3];
rz(-1.1223464) q[3];
sx q[3];
rz(-0.34235172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0943429) q[2];
sx q[2];
rz(-0.80545682) q[2];
sx q[2];
rz(-0.94364014) q[2];
rz(2.0761679) q[3];
sx q[3];
rz(-1.7908955) q[3];
sx q[3];
rz(2.0431199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074325398) q[0];
sx q[0];
rz(-1.704498) q[0];
sx q[0];
rz(-0.0083010439) q[0];
rz(-0.51757327) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(1.5932721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2733214) q[0];
sx q[0];
rz(-1.5948606) q[0];
sx q[0];
rz(-1.2198183) q[0];
x q[1];
rz(-0.50524925) q[2];
sx q[2];
rz(-1.8191686) q[2];
sx q[2];
rz(2.5587683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1344117) q[1];
sx q[1];
rz(-0.500965) q[1];
sx q[1];
rz(0.8207313) q[1];
x q[2];
rz(-2.1861784) q[3];
sx q[3];
rz(-1.9166167) q[3];
sx q[3];
rz(2.7503848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3090618) q[2];
sx q[2];
rz(-1.6522202) q[2];
sx q[2];
rz(1.8050516) q[2];
rz(3.0858223) q[3];
sx q[3];
rz(-0.99168188) q[3];
sx q[3];
rz(-2.706004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5830773) q[0];
sx q[0];
rz(-0.4774839) q[0];
sx q[0];
rz(0.68786311) q[0];
rz(1.4631924) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(-2.8942143) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32175049) q[0];
sx q[0];
rz(-2.8341056) q[0];
sx q[0];
rz(-0.64096682) q[0];
x q[1];
rz(-2.3937463) q[2];
sx q[2];
rz(-2.3299542) q[2];
sx q[2];
rz(-2.5876665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0717031) q[1];
sx q[1];
rz(-2.1462198) q[1];
sx q[1];
rz(-0.50168049) q[1];
rz(1.803627) q[3];
sx q[3];
rz(-1.4484753) q[3];
sx q[3];
rz(-1.8450774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4443724) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(0.42416254) q[2];
rz(2.0992725) q[3];
sx q[3];
rz(-2.3764231) q[3];
sx q[3];
rz(2.585129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0031849) q[0];
sx q[0];
rz(-2.1557032) q[0];
sx q[0];
rz(2.5307122) q[0];
rz(-2.8914087) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(2.6649323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7358897) q[0];
sx q[0];
rz(-2.6578356) q[0];
sx q[0];
rz(-1.4866845) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0518529) q[2];
sx q[2];
rz(-1.753627) q[2];
sx q[2];
rz(0.51328608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5077153) q[1];
sx q[1];
rz(-2.4832279) q[1];
sx q[1];
rz(-2.4509199) q[1];
x q[2];
rz(0.96673217) q[3];
sx q[3];
rz(-2.4632404) q[3];
sx q[3];
rz(-2.683992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99159795) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-0.67982137) q[2];
rz(2.6796135) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(1.6405039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221472) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(-2.9875901) q[0];
rz(2.9601861) q[1];
sx q[1];
rz(-1.0702952) q[1];
sx q[1];
rz(-3.0214686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87163602) q[0];
sx q[0];
rz(-1.1678732) q[0];
sx q[0];
rz(-0.064489207) q[0];
rz(-pi) q[1];
rz(0.93379069) q[2];
sx q[2];
rz(-2.8178038) q[2];
sx q[2];
rz(2.0227014) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57615818) q[1];
sx q[1];
rz(-1.6685969) q[1];
sx q[1];
rz(1.3937772) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8974344) q[3];
sx q[3];
rz(-2.6006581) q[3];
sx q[3];
rz(1.8997826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4274365) q[2];
sx q[2];
rz(-1.7566046) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(0.51042026) q[3];
sx q[3];
rz(-2.3000058) q[3];
sx q[3];
rz(-0.69971219) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0521312) q[0];
sx q[0];
rz(-2.1928146) q[0];
sx q[0];
rz(-2.7556038) q[0];
rz(-1.5929068) q[1];
sx q[1];
rz(-2.6451151) q[1];
sx q[1];
rz(-1.5923502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44819866) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(1.4017808) q[0];
rz(2.4108294) q[2];
sx q[2];
rz(-2.363702) q[2];
sx q[2];
rz(-0.90532263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85268116) q[1];
sx q[1];
rz(-1.4687454) q[1];
sx q[1];
rz(-0.026983807) q[1];
rz(-pi) q[2];
rz(0.066927197) q[3];
sx q[3];
rz(-1.7239583) q[3];
sx q[3];
rz(-2.7471971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7768895) q[2];
sx q[2];
rz(-0.93901712) q[2];
sx q[2];
rz(1.6949863) q[2];
rz(-2.1052836) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(0.026329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885289) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(-0.22008315) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(-1.2575116) q[2];
sx q[2];
rz(-0.85713119) q[2];
sx q[2];
rz(-1.5766889) q[2];
rz(0.50095273) q[3];
sx q[3];
rz(-1.2798449) q[3];
sx q[3];
rz(2.4191054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
