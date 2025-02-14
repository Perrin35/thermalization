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
rz(0.96639291) q[0];
sx q[0];
rz(-1.7557431) q[0];
sx q[0];
rz(-0.59544271) q[0];
rz(1.2915986) q[1];
sx q[1];
rz(-2.5505677) q[1];
sx q[1];
rz(-0.12330595) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46245136) q[0];
sx q[0];
rz(-1.9872905) q[0];
sx q[0];
rz(0.97521675) q[0];
rz(0.79084556) q[2];
sx q[2];
rz(-0.95480761) q[2];
sx q[2];
rz(0.14033422) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1327269) q[1];
sx q[1];
rz(-1.4767329) q[1];
sx q[1];
rz(-1.3558741) q[1];
rz(-pi) q[2];
rz(1.7233347) q[3];
sx q[3];
rz(-1.2087052) q[3];
sx q[3];
rz(0.83521508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3143602) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(-3.0296791) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.1575969) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258485) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-0.14990212) q[0];
rz(2.8556178) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(-1.1289977) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1492831) q[0];
sx q[0];
rz(-0.2899) q[0];
sx q[0];
rz(-1.4647746) q[0];
rz(-pi) q[1];
rz(1.6304134) q[2];
sx q[2];
rz(-1.5853512) q[2];
sx q[2];
rz(-2.3090135) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2098238) q[1];
sx q[1];
rz(-1.9289513) q[1];
sx q[1];
rz(-1.7478554) q[1];
rz(-pi) q[2];
rz(-0.942248) q[3];
sx q[3];
rz(-1.0744922) q[3];
sx q[3];
rz(1.4013578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67843208) q[2];
sx q[2];
rz(-0.9223991) q[2];
sx q[2];
rz(-1.9909667) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(-0.020603389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(-2.7599957) q[0];
rz(-1.8474139) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(0.19827422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0509153) q[0];
sx q[0];
rz(-0.28813513) q[0];
sx q[0];
rz(-2.227562) q[0];
rz(-pi) q[1];
rz(-0.080520544) q[2];
sx q[2];
rz(-0.90290194) q[2];
sx q[2];
rz(0.59216796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5845909) q[1];
sx q[1];
rz(-1.7832568) q[1];
sx q[1];
rz(0.23834385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4573077) q[3];
sx q[3];
rz(-0.17593613) q[3];
sx q[3];
rz(0.709155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76564378) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(0.51885968) q[2];
rz(1.8910003) q[3];
sx q[3];
rz(-2.0542681) q[3];
sx q[3];
rz(-0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4984703) q[0];
sx q[0];
rz(-2.9415218) q[0];
sx q[0];
rz(0.44922391) q[0];
rz(1.6953702) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(-0.62072388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0064471) q[0];
sx q[0];
rz(-0.9165316) q[0];
sx q[0];
rz(-1.7530789) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5081257) q[2];
sx q[2];
rz(-1.4711498) q[2];
sx q[2];
rz(1.5657305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52445946) q[1];
sx q[1];
rz(-1.7726328) q[1];
sx q[1];
rz(2.733888) q[1];
rz(-3.0884864) q[3];
sx q[3];
rz(-1.0328919) q[3];
sx q[3];
rz(-0.90467473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(-0.16109666) q[2];
rz(-0.39786878) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(-0.14959344) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38254) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-0.6404666) q[1];
sx q[1];
rz(1.2535198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81957626) q[0];
sx q[0];
rz(-1.6068234) q[0];
sx q[0];
rz(3.1234804) q[0];
x q[1];
rz(1.8158556) q[2];
sx q[2];
rz(-0.7759717) q[2];
sx q[2];
rz(1.4910335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1977928) q[1];
sx q[1];
rz(-1.6159004) q[1];
sx q[1];
rz(0.79774858) q[1];
x q[2];
rz(-0.34714602) q[3];
sx q[3];
rz(-1.5724564) q[3];
sx q[3];
rz(2.0059137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(1.8750635) q[2];
rz(0.62147102) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(-2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(0.025064502) q[0];
rz(1.2114245) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(2.8866344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0421521) q[0];
sx q[0];
rz(-1.6201311) q[0];
sx q[0];
rz(1.5758324) q[0];
x q[1];
rz(2.4676085) q[2];
sx q[2];
rz(-2.4366798) q[2];
sx q[2];
rz(2.2121536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7656169) q[1];
sx q[1];
rz(-2.097953) q[1];
sx q[1];
rz(-0.65816452) q[1];
x q[2];
rz(-0.21277748) q[3];
sx q[3];
rz(-1.0085953) q[3];
sx q[3];
rz(-0.24383814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-1.0422372) q[2];
sx q[2];
rz(-0.45905217) q[2];
rz(-1.6445271) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(3.0204401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727305) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(-0.91880265) q[1];
sx q[1];
rz(-2.0252392) q[1];
sx q[1];
rz(1.2109717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54706193) q[0];
sx q[0];
rz(-1.8939202) q[0];
sx q[0];
rz(1.7842152) q[0];
x q[1];
rz(2.2077256) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(-3.0702555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0414435) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(-0.99975296) q[1];
rz(-pi) q[2];
rz(-1.2065229) q[3];
sx q[3];
rz(-1.9134054) q[3];
sx q[3];
rz(-0.87876696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81898895) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(-0.71713478) q[2];
rz(2.1142193) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-0.43924847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1502007) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(3.126934) q[0];
rz(-2.390059) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(-1.6709447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6122307) q[0];
sx q[0];
rz(-0.23461537) q[0];
sx q[0];
rz(0.63572065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8445817) q[2];
sx q[2];
rz(-2.8737465) q[2];
sx q[2];
rz(1.0849407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60622207) q[1];
sx q[1];
rz(-0.41946966) q[1];
sx q[1];
rz(1.8505881) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1378391) q[3];
sx q[3];
rz(-1.9158746) q[3];
sx q[3];
rz(1.0964637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26059255) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(0.27349681) q[2];
rz(1.986844) q[3];
sx q[3];
rz(-0.65360779) q[3];
sx q[3];
rz(-0.65034136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1028033) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(1.0754841) q[0];
rz(0.12818809) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(0.34559616) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81886473) q[0];
sx q[0];
rz(-2.4562307) q[0];
sx q[0];
rz(-1.8249874) q[0];
rz(1.0900201) q[2];
sx q[2];
rz(-0.62453237) q[2];
sx q[2];
rz(-0.20779729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84050485) q[1];
sx q[1];
rz(-1.6337758) q[1];
sx q[1];
rz(-2.5999864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4050499) q[3];
sx q[3];
rz(-0.9791383) q[3];
sx q[3];
rz(1.4480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0868316) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(-1.5209939) q[2];
rz(-1.7704376) q[3];
sx q[3];
rz(-2.3833279) q[3];
sx q[3];
rz(-0.41485205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.82776752) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-0.23571043) q[0];
rz(2.56855) q[1];
sx q[1];
rz(-2.133281) q[1];
sx q[1];
rz(1.3689573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9873857) q[0];
sx q[0];
rz(-1.1416178) q[0];
sx q[0];
rz(2.3230419) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4659285) q[2];
sx q[2];
rz(-1.079529) q[2];
sx q[2];
rz(-1.5671687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6073709) q[1];
sx q[1];
rz(-1.0080308) q[1];
sx q[1];
rz(-2.9396179) q[1];
rz(-pi) q[2];
rz(-0.36054109) q[3];
sx q[3];
rz(-1.6586132) q[3];
sx q[3];
rz(-0.58780625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3706751) q[2];
sx q[2];
rz(-1.8149899) q[2];
sx q[2];
rz(3.0885922) q[2];
rz(0.75183374) q[3];
sx q[3];
rz(-2.1078608) q[3];
sx q[3];
rz(-0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.62405217) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(1.7908295) q[1];
sx q[1];
rz(-2.7443934) q[1];
sx q[1];
rz(2.9846356) q[1];
rz(-0.69722947) q[2];
sx q[2];
rz(-1.4705428) q[2];
sx q[2];
rz(-1.85208) q[2];
rz(-0.80812576) q[3];
sx q[3];
rz(-1.1302409) q[3];
sx q[3];
rz(-1.7779868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
