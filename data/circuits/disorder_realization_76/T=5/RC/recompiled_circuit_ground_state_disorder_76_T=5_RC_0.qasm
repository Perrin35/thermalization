OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4467093) q[0];
sx q[0];
rz(-2.0599685) q[0];
sx q[0];
rz(-2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(1.4256328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6524871) q[0];
sx q[0];
rz(-1.8579626) q[0];
sx q[0];
rz(1.0190359) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6447634) q[2];
sx q[2];
rz(-1.8112) q[2];
sx q[2];
rz(-2.5508326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.718218) q[1];
sx q[1];
rz(-1.4932695) q[1];
sx q[1];
rz(-1.7192057) q[1];
rz(-pi) q[2];
rz(-2.2896122) q[3];
sx q[3];
rz(-1.6036828) q[3];
sx q[3];
rz(-2.6926958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2322959) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(-1.9624814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772684) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(2.6569195) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.1601123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697333) q[0];
sx q[0];
rz(-1.5892649) q[0];
sx q[0];
rz(2.8375677) q[0];
rz(-2.2641029) q[2];
sx q[2];
rz(-0.98131591) q[2];
sx q[2];
rz(1.9716919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96148032) q[1];
sx q[1];
rz(-2.4870076) q[1];
sx q[1];
rz(-2.727319) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6187702) q[3];
sx q[3];
rz(-1.8588052) q[3];
sx q[3];
rz(-3.0436181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(0.45787946) q[2];
rz(-2.0186021) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(2.5751953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(-0.031524468) q[0];
rz(2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.9256176) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77338868) q[0];
sx q[0];
rz(-0.88951123) q[0];
sx q[0];
rz(0.70685203) q[0];
x q[1];
rz(-2.6570733) q[2];
sx q[2];
rz(-1.7555825) q[2];
sx q[2];
rz(-0.69714025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2906704) q[1];
sx q[1];
rz(-1.8462204) q[1];
sx q[1];
rz(-1.2864134) q[1];
rz(-pi) q[2];
rz(-2.2103045) q[3];
sx q[3];
rz(-2.1289325) q[3];
sx q[3];
rz(-1.0071514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1027801) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(-0.83596027) q[2];
rz(1.7838259) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(-0.74762216) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589856) q[0];
sx q[0];
rz(-1.405412) q[0];
sx q[0];
rz(-0.080168515) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(1.3287883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565312) q[0];
sx q[0];
rz(-2.079112) q[0];
sx q[0];
rz(-1.9446745) q[0];
rz(-pi) q[1];
rz(-2.8187739) q[2];
sx q[2];
rz(-1.7344966) q[2];
sx q[2];
rz(1.9893008) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7419974) q[1];
sx q[1];
rz(-2.6873702) q[1];
sx q[1];
rz(1.2356367) q[1];
x q[2];
rz(2.8058047) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(1.646281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(2.4713016) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(-1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9012673) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(-1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-2.7010837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021316) q[0];
sx q[0];
rz(-1.4533784) q[0];
sx q[0];
rz(1.1068547) q[0];
rz(-0.20655234) q[2];
sx q[2];
rz(-1.9777386) q[2];
sx q[2];
rz(-0.47688866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9677093) q[1];
sx q[1];
rz(-1.5313799) q[1];
sx q[1];
rz(2.8110678) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83310762) q[3];
sx q[3];
rz(-0.78189497) q[3];
sx q[3];
rz(1.1469008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0195007) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(-2.7276373) q[2];
rz(-1.7533938) q[3];
sx q[3];
rz(-1.5475169) q[3];
sx q[3];
rz(-0.12510124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6269161) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(0.93389121) q[0];
rz(2.4712708) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(1.3353039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437296) q[0];
sx q[0];
rz(-1.4708733) q[0];
sx q[0];
rz(-2.1735682) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.983903) q[2];
sx q[2];
rz(-1.5458428) q[2];
sx q[2];
rz(1.4098997) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.858243) q[1];
sx q[1];
rz(-1.5482386) q[1];
sx q[1];
rz(1.4580887) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1298529) q[3];
sx q[3];
rz(-0.63705963) q[3];
sx q[3];
rz(0.54367346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(-1.7719841) q[3];
sx q[3];
rz(-1.2415875) q[3];
sx q[3];
rz(-0.52186596) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(-0.15368803) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.2362365) q[1];
sx q[1];
rz(2.1930146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0517387) q[0];
sx q[0];
rz(-1.3696545) q[0];
sx q[0];
rz(-0.62369831) q[0];
rz(-pi) q[1];
rz(-0.3237299) q[2];
sx q[2];
rz(-2.1846131) q[2];
sx q[2];
rz(0.15440369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7708009) q[1];
sx q[1];
rz(-1.3389412) q[1];
sx q[1];
rz(-1.8370017) q[1];
rz(-1.7465215) q[3];
sx q[3];
rz(-2.4843289) q[3];
sx q[3];
rz(0.45096179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1367246) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(-3.1401805) q[2];
rz(-3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(0.96127659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.75981265) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(2.5835999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787631) q[0];
sx q[0];
rz(-0.42233322) q[0];
sx q[0];
rz(-0.93393737) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9079014) q[2];
sx q[2];
rz(-0.56392852) q[2];
sx q[2];
rz(0.83966161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26436603) q[1];
sx q[1];
rz(-2.5747882) q[1];
sx q[1];
rz(1.0344489) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6648983) q[3];
sx q[3];
rz(-0.27827874) q[3];
sx q[3];
rz(-0.32839963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(0.26926678) q[2];
rz(-1.1632129) q[3];
sx q[3];
rz(-0.60267699) q[3];
sx q[3];
rz(-2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.4686541) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-0.55666298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024324) q[0];
sx q[0];
rz(-0.35050979) q[0];
sx q[0];
rz(-1.9950161) q[0];
x q[1];
rz(0.78293856) q[2];
sx q[2];
rz(-0.73854827) q[2];
sx q[2];
rz(2.1434181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4924188) q[1];
sx q[1];
rz(-0.88090501) q[1];
sx q[1];
rz(-1.0997195) q[1];
rz(1.6708371) q[3];
sx q[3];
rz(-1.8042068) q[3];
sx q[3];
rz(1.4518713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2399981) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(1.979801) q[2];
rz(0.53317541) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(-0.1575135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.6530782) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(-2.5467806) q[0];
rz(-0.57890233) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(-2.4822809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4730277) q[0];
sx q[0];
rz(-2.1921232) q[0];
sx q[0];
rz(-1.4061808) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5388902) q[2];
sx q[2];
rz(-1.1991767) q[2];
sx q[2];
rz(-2.0211257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.59567) q[1];
sx q[1];
rz(-2.8575142) q[1];
sx q[1];
rz(2.0998663) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80307428) q[3];
sx q[3];
rz(-0.37130203) q[3];
sx q[3];
rz(0.40483958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.351563) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(3.0885546) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(-0.0742577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.2863083) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(3.042649) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(-3.121539) q[2];
sx q[2];
rz(-0.40553025) q[2];
sx q[2];
rz(2.5237906) q[2];
rz(2.4002038) q[3];
sx q[3];
rz(-1.4546053) q[3];
sx q[3];
rz(-1.4447053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
