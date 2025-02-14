OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.111423) q[0];
sx q[0];
rz(-1.3601967) q[0];
sx q[0];
rz(-2.7786322) q[0];
rz(1.9650004) q[1];
sx q[1];
rz(-1.3624374) q[1];
sx q[1];
rz(-1.6972313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43508726) q[0];
sx q[0];
rz(-1.6643063) q[0];
sx q[0];
rz(2.963752) q[0];
rz(-pi) q[1];
rz(-1.6955201) q[2];
sx q[2];
rz(-1.8565923) q[2];
sx q[2];
rz(0.054626183) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.23420424) q[1];
sx q[1];
rz(-1.475913) q[1];
sx q[1];
rz(2.0070875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6300171) q[3];
sx q[3];
rz(-2.539336) q[3];
sx q[3];
rz(1.2685552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1732257) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(3.0464029) q[2];
rz(1.7732874) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(-0.29909721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.478941) q[0];
sx q[0];
rz(-0.79942411) q[0];
sx q[0];
rz(2.4484402) q[0];
rz(-2.8043546) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(0.79280028) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5186441) q[0];
sx q[0];
rz(-2.3272133) q[0];
sx q[0];
rz(-1.9163535) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9117965) q[2];
sx q[2];
rz(-2.0760798) q[2];
sx q[2];
rz(-1.2838703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2094272) q[1];
sx q[1];
rz(-1.2840198) q[1];
sx q[1];
rz(0.78358235) q[1];
rz(-pi) q[2];
rz(-0.81962193) q[3];
sx q[3];
rz(-1.9385129) q[3];
sx q[3];
rz(-1.6727656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99337951) q[2];
sx q[2];
rz(-2.9714163) q[2];
sx q[2];
rz(-1.4519838) q[2];
rz(-0.65886894) q[3];
sx q[3];
rz(-1.6784607) q[3];
sx q[3];
rz(0.39670593) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0186998) q[0];
sx q[0];
rz(-1.0188228) q[0];
sx q[0];
rz(3.0320211) q[0];
rz(-2.295629) q[1];
sx q[1];
rz(-0.94596091) q[1];
sx q[1];
rz(-0.74839655) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0349265) q[0];
sx q[0];
rz(-2.1252118) q[0];
sx q[0];
rz(-3.0361452) q[0];
x q[1];
rz(-0.72785925) q[2];
sx q[2];
rz(-2.4185491) q[2];
sx q[2];
rz(-2.7287908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8037532) q[1];
sx q[1];
rz(-0.0633792) q[1];
sx q[1];
rz(1.4638639) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.844039) q[3];
sx q[3];
rz(-0.58876172) q[3];
sx q[3];
rz(-1.2514284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1154068) q[2];
sx q[2];
rz(-1.7188965) q[2];
sx q[2];
rz(2.4288948) q[2];
rz(1.4767492) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(-0.084499806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0067714) q[0];
sx q[0];
rz(-1.1399784) q[0];
sx q[0];
rz(-0.086932927) q[0];
rz(0.25761071) q[1];
sx q[1];
rz(-1.0898217) q[1];
sx q[1];
rz(1.2923406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564926) q[0];
sx q[0];
rz(-0.79946639) q[0];
sx q[0];
rz(1.234458) q[0];
x q[1];
rz(-1.7316225) q[2];
sx q[2];
rz(-1.8984744) q[2];
sx q[2];
rz(0.92152061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.41117696) q[1];
sx q[1];
rz(-2.431793) q[1];
sx q[1];
rz(2.1201154) q[1];
rz(2.7601808) q[3];
sx q[3];
rz(-2.7578045) q[3];
sx q[3];
rz(-0.82510766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8452235) q[2];
sx q[2];
rz(-1.951428) q[2];
sx q[2];
rz(-0.84687084) q[2];
rz(1.3891247) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(-2.5510767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95319372) q[0];
sx q[0];
rz(-2.0162835) q[0];
sx q[0];
rz(-0.97287792) q[0];
rz(2.1994622) q[1];
sx q[1];
rz(-0.82754389) q[1];
sx q[1];
rz(-0.00076248893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1050232) q[0];
sx q[0];
rz(-2.3876086) q[0];
sx q[0];
rz(-0.79490173) q[0];
x q[1];
rz(2.7050651) q[2];
sx q[2];
rz(-0.68219705) q[2];
sx q[2];
rz(1.4774708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.768133) q[1];
sx q[1];
rz(-0.64300696) q[1];
sx q[1];
rz(2.0483584) q[1];
x q[2];
rz(0.11526626) q[3];
sx q[3];
rz(-2.5433833) q[3];
sx q[3];
rz(-2.0245352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7611277) q[2];
sx q[2];
rz(-1.7385812) q[2];
sx q[2];
rz(-2.7565956) q[2];
rz(-0.72832406) q[3];
sx q[3];
rz(-0.60898048) q[3];
sx q[3];
rz(-2.716632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4094792) q[0];
sx q[0];
rz(-2.3347169) q[0];
sx q[0];
rz(-2.683486) q[0];
rz(-1.7181646) q[1];
sx q[1];
rz(-0.79045311) q[1];
sx q[1];
rz(3.0480393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76339611) q[0];
sx q[0];
rz(-0.79471055) q[0];
sx q[0];
rz(2.2346157) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92624591) q[2];
sx q[2];
rz(-2.8684596) q[2];
sx q[2];
rz(-0.63970837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72274745) q[1];
sx q[1];
rz(-1.0422575) q[1];
sx q[1];
rz(1.8209807) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5922067) q[3];
sx q[3];
rz(-1.011102) q[3];
sx q[3];
rz(-2.7407854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88390049) q[2];
sx q[2];
rz(-0.42282405) q[2];
sx q[2];
rz(2.532393) q[2];
rz(0.60112634) q[3];
sx q[3];
rz(-1.5355891) q[3];
sx q[3];
rz(0.48457423) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5179829) q[0];
sx q[0];
rz(-1.6684775) q[0];
sx q[0];
rz(1.8286888) q[0];
rz(3.0598705) q[1];
sx q[1];
rz(-1.7535968) q[1];
sx q[1];
rz(2.0770226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3400187) q[0];
sx q[0];
rz(-0.94807762) q[0];
sx q[0];
rz(2.4452371) q[0];
rz(2.1871952) q[2];
sx q[2];
rz(-1.5845044) q[2];
sx q[2];
rz(-0.15132667) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7408167) q[1];
sx q[1];
rz(-2.1451752) q[1];
sx q[1];
rz(-3.0568491) q[1];
rz(0.40708159) q[3];
sx q[3];
rz(-1.0069869) q[3];
sx q[3];
rz(-2.8460007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99819034) q[2];
sx q[2];
rz(-2.7394131) q[2];
sx q[2];
rz(1.1581988) q[2];
rz(-2.4573333) q[3];
sx q[3];
rz(-1.8161769) q[3];
sx q[3];
rz(-1.6392596) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96712464) q[0];
sx q[0];
rz(-0.10395771) q[0];
sx q[0];
rz(0.52484584) q[0];
rz(1.9606494) q[1];
sx q[1];
rz(-2.3044105) q[1];
sx q[1];
rz(0.7276853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8129163) q[0];
sx q[0];
rz(-0.70309454) q[0];
sx q[0];
rz(1.9090184) q[0];
rz(-pi) q[1];
rz(1.069017) q[2];
sx q[2];
rz(-2.2282218) q[2];
sx q[2];
rz(1.1595961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3062417) q[1];
sx q[1];
rz(-2.0064276) q[1];
sx q[1];
rz(-1.8620051) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0342909) q[3];
sx q[3];
rz(-0.91917097) q[3];
sx q[3];
rz(-0.95355031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9119447) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(-1.7952807) q[2];
rz(0.66156578) q[3];
sx q[3];
rz(-1.1450359) q[3];
sx q[3];
rz(-0.0865817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10053703) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(2.2030785) q[0];
rz(1.8427294) q[1];
sx q[1];
rz(-0.88911903) q[1];
sx q[1];
rz(-3.0079957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.035048) q[0];
sx q[0];
rz(-1.4770916) q[0];
sx q[0];
rz(0.32982112) q[0];
rz(-2.5659525) q[2];
sx q[2];
rz(-1.8274022) q[2];
sx q[2];
rz(2.9413066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6386968) q[1];
sx q[1];
rz(-1.2112482) q[1];
sx q[1];
rz(-1.6537958) q[1];
x q[2];
rz(1.2032005) q[3];
sx q[3];
rz(-2.1798686) q[3];
sx q[3];
rz(-1.4673135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9962697) q[2];
sx q[2];
rz(-2.1349553) q[2];
sx q[2];
rz(0.88627306) q[2];
rz(-0.26850548) q[3];
sx q[3];
rz(-2.0566745) q[3];
sx q[3];
rz(1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71235424) q[0];
sx q[0];
rz(-0.3952643) q[0];
sx q[0];
rz(0.23671737) q[0];
rz(0.16353823) q[1];
sx q[1];
rz(-1.5312559) q[1];
sx q[1];
rz(-2.9530361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9360787) q[0];
sx q[0];
rz(-0.88919836) q[0];
sx q[0];
rz(-2.71191) q[0];
x q[1];
rz(-0.46444664) q[2];
sx q[2];
rz(-2.6147343) q[2];
sx q[2];
rz(1.3642481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8324667) q[1];
sx q[1];
rz(-1.6285076) q[1];
sx q[1];
rz(-2.0821435) q[1];
rz(-pi) q[2];
rz(0.43980912) q[3];
sx q[3];
rz(-1.4521003) q[3];
sx q[3];
rz(-2.4921093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7207555) q[2];
sx q[2];
rz(-1.5034224) q[2];
sx q[2];
rz(2.3756964) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(-2.6115821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6388549) q[0];
sx q[0];
rz(-1.9866332) q[0];
sx q[0];
rz(-1.6992983) q[0];
rz(2.8522708) q[1];
sx q[1];
rz(-1.0565636) q[1];
sx q[1];
rz(0.4531959) q[1];
rz(2.6027457) q[2];
sx q[2];
rz(-1.4543903) q[2];
sx q[2];
rz(-1.7001154) q[2];
rz(3.0719793) q[3];
sx q[3];
rz(-1.1021464) q[3];
sx q[3];
rz(-2.7819421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
