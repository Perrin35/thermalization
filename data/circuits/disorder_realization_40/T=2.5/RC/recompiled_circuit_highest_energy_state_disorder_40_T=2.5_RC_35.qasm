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
rz(1.9873729) q[0];
sx q[0];
rz(4.6648751) q[0];
sx q[0];
rz(9.2821791) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(3.7945336) q[1];
sx q[1];
rz(8.3795587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088974911) q[0];
sx q[0];
rz(-1.0427022) q[0];
sx q[0];
rz(-2.108197) q[0];
rz(-pi) q[1];
rz(-3.1047761) q[2];
sx q[2];
rz(-2.2082001) q[2];
sx q[2];
rz(0.76257818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(-1.0008903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0020489) q[3];
sx q[3];
rz(-1.6709329) q[3];
sx q[3];
rz(1.9619693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1656701) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(2.9284076) q[2];
rz(-2.2255157) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(0.38755125) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857472) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(0.47072738) q[0];
rz(-1.1031021) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(-2.5618166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743917) q[0];
sx q[0];
rz(-1.5109332) q[0];
sx q[0];
rz(-0.52459985) q[0];
rz(-pi) q[1];
rz(0.51885278) q[2];
sx q[2];
rz(-2.5343072) q[2];
sx q[2];
rz(0.80087207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3892711) q[1];
sx q[1];
rz(-1.0022517) q[1];
sx q[1];
rz(2.6804377) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2006496) q[3];
sx q[3];
rz(-1.3251628) q[3];
sx q[3];
rz(0.64526886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(-0.095005438) q[2];
rz(-2.5659918) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.59548241) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(0.27221671) q[0];
rz(0.1046003) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(1.4629755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7143216) q[0];
sx q[0];
rz(-2.2939689) q[0];
sx q[0];
rz(-3.0747674) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5065423) q[2];
sx q[2];
rz(-2.3550866) q[2];
sx q[2];
rz(0.74555874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.453664) q[1];
sx q[1];
rz(-1.8519823) q[1];
sx q[1];
rz(2.4635876) q[1];
rz(-pi) q[2];
rz(-1.6312863) q[3];
sx q[3];
rz(-1.9391427) q[3];
sx q[3];
rz(0.18243901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0649123) q[2];
sx q[2];
rz(-1.400759) q[2];
sx q[2];
rz(-2.7317375) q[2];
rz(-0.05154933) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16875295) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(-2.4421316) q[0];
rz(-0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(1.0994937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52877766) q[0];
sx q[0];
rz(-1.8474425) q[0];
sx q[0];
rz(-2.1601281) q[0];
rz(-pi) q[1];
rz(-1.0223939) q[2];
sx q[2];
rz(-1.7779997) q[2];
sx q[2];
rz(-2.0275379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.857497) q[1];
sx q[1];
rz(-1.2624143) q[1];
sx q[1];
rz(0.33607884) q[1];
rz(-pi) q[2];
rz(-0.52249281) q[3];
sx q[3];
rz(-0.71194369) q[3];
sx q[3];
rz(-0.57761907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76116556) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(-0.94044828) q[2];
rz(-2.579651) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(-0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024427323) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(1.0166136) q[0];
rz(0.92574614) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0177855) q[0];
sx q[0];
rz(-1.9636256) q[0];
sx q[0];
rz(0.17946243) q[0];
rz(-0.87575298) q[2];
sx q[2];
rz(-0.31538439) q[2];
sx q[2];
rz(-1.8502667) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13097732) q[1];
sx q[1];
rz(-1.0189799) q[1];
sx q[1];
rz(0.24893649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29239817) q[3];
sx q[3];
rz(-0.83857036) q[3];
sx q[3];
rz(-0.75493073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(3.0184271) q[2];
rz(2.4672616) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(-0.82790747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8318361) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(-0.24093974) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(0.13490881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19770007) q[0];
sx q[0];
rz(-2.0257844) q[0];
sx q[0];
rz(0.83857341) q[0];
rz(-pi) q[1];
rz(0.1730072) q[2];
sx q[2];
rz(-2.5998625) q[2];
sx q[2];
rz(2.2587905) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.125306) q[1];
sx q[1];
rz(-1.8452438) q[1];
sx q[1];
rz(-2.3322565) q[1];
rz(-0.44996275) q[3];
sx q[3];
rz(-1.5057766) q[3];
sx q[3];
rz(-1.3715594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2065108) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(-1.8180465) q[2];
rz(2.8694618) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(-2.4830699) q[0];
rz(1.5918484) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(0.50643593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69909912) q[0];
sx q[0];
rz(-2.7662219) q[0];
sx q[0];
rz(0.18113776) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6295048) q[2];
sx q[2];
rz(-0.48248267) q[2];
sx q[2];
rz(0.62472945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68489198) q[1];
sx q[1];
rz(-1.4948256) q[1];
sx q[1];
rz(3.0481245) q[1];
rz(0.063252216) q[3];
sx q[3];
rz(-0.58025415) q[3];
sx q[3];
rz(-1.3903416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3353614) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.7757724) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34614554) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(-0.090713352) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(1.027164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5266383) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(2.5102708) q[0];
rz(-pi) q[1];
x q[1];
rz(0.022826957) q[2];
sx q[2];
rz(-1.5787303) q[2];
sx q[2];
rz(1.4102907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5324983) q[1];
sx q[1];
rz(-0.23819345) q[1];
sx q[1];
rz(1.1042117) q[1];
x q[2];
rz(1.5040565) q[3];
sx q[3];
rz(-0.68664521) q[3];
sx q[3];
rz(-2.8444169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(0.99009222) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(0.038343553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7541499) q[0];
sx q[0];
rz(-0.70873547) q[0];
sx q[0];
rz(-0.55737525) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(0.16709669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1633269) q[0];
sx q[0];
rz(-1.3045132) q[0];
sx q[0];
rz(2.8243218) q[0];
rz(-2.4403247) q[2];
sx q[2];
rz(-1.7116364) q[2];
sx q[2];
rz(0.62447157) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5846338) q[1];
sx q[1];
rz(-1.2642197) q[1];
sx q[1];
rz(-2.4584998) q[1];
rz(-pi) q[2];
rz(-1.6623508) q[3];
sx q[3];
rz(-2.4560438) q[3];
sx q[3];
rz(0.26622546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(2.7131405) q[2];
rz(-2.4109449) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(-2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927602) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(1.1195419) q[0];
rz(-0.66850942) q[1];
sx q[1];
rz(-0.57066494) q[1];
sx q[1];
rz(-2.8817435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90454532) q[0];
sx q[0];
rz(-2.0656842) q[0];
sx q[0];
rz(-0.2965692) q[0];
rz(-pi) q[1];
rz(1.1041743) q[2];
sx q[2];
rz(-1.9167308) q[2];
sx q[2];
rz(-0.26901252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6653241) q[1];
sx q[1];
rz(-0.55108738) q[1];
sx q[1];
rz(2.4399806) q[1];
rz(-pi) q[2];
rz(2.9838461) q[3];
sx q[3];
rz(-2.254527) q[3];
sx q[3];
rz(0.58303787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4708289) q[2];
sx q[2];
rz(-1.2329817) q[2];
sx q[2];
rz(-2.0551576) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-2.2958675) q[3];
sx q[3];
rz(0.72211784) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54871854) q[0];
sx q[0];
rz(-1.6328136) q[0];
sx q[0];
rz(1.6528224) q[0];
rz(-1.2552352) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(1.3484501) q[2];
sx q[2];
rz(-2.7334474) q[2];
sx q[2];
rz(1.853142) q[2];
rz(1.6519573) q[3];
sx q[3];
rz(-2.4774144) q[3];
sx q[3];
rz(-1.7437205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
