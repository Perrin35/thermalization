OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(2.6775223) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(-2.2201559) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48972763) q[2];
sx q[2];
rz(-1.5416607) q[2];
sx q[2];
rz(2.4907128) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97548188) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(2.6544177) q[1];
rz(-0.10847096) q[3];
sx q[3];
rz(-0.15158187) q[3];
sx q[3];
rz(0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-0.79663509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684727) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(-0.40533439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7700427) q[2];
sx q[2];
rz(-2.1363878) q[2];
sx q[2];
rz(-2.9233962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4371725) q[1];
sx q[1];
rz(-2.17802) q[1];
sx q[1];
rz(1.5551268) q[1];
x q[2];
rz(0.96062406) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(2.744439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1709135) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(-0.12165102) q[0];
x q[1];
rz(0.59962745) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(-1.8011013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2726047) q[1];
sx q[1];
rz(-1.4613978) q[1];
sx q[1];
rz(-1.617856) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(0.56328303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3455428) q[0];
sx q[0];
rz(-0.69957083) q[0];
sx q[0];
rz(2.0027341) q[0];
rz(1.9070542) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(-1.3627571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(2.8698189) q[1];
x q[2];
rz(1.8654278) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(-2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974385) q[0];
sx q[0];
rz(-2.4442406) q[0];
sx q[0];
rz(-2.5028412) q[0];
rz(-pi) q[1];
rz(-2.7439793) q[2];
sx q[2];
rz(-2.0587454) q[2];
sx q[2];
rz(0.60148009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3106766) q[1];
sx q[1];
rz(-1.1569996) q[1];
sx q[1];
rz(0.78891854) q[1];
rz(-0.01708548) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(-2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(-0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.8355339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40697843) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(2.2218496) q[0];
x q[1];
rz(1.3681709) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(-2.1855598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8514511) q[1];
sx q[1];
rz(-0.72039225) q[1];
sx q[1];
rz(3.1300701) q[1];
rz(-pi) q[2];
rz(1.7197051) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(1.0466928) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55506599) q[0];
sx q[0];
rz(-1.5952206) q[0];
sx q[0];
rz(1.3506372) q[0];
x q[1];
rz(-2.3796758) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(-2.6921536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5187832) q[1];
sx q[1];
rz(-0.86474027) q[1];
sx q[1];
rz(-2.519636) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8700637) q[3];
sx q[3];
rz(-0.93942681) q[3];
sx q[3];
rz(2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053514078) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(1.0900351) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8939395) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(2.0695956) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1689414) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0197969) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(1.4956054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-0.91167489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.002279) q[0];
sx q[0];
rz(-0.3542491) q[0];
sx q[0];
rz(-0.87689633) q[0];
rz(-1.6976835) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(-1.1586231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4166491) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(-1.7979421) q[1];
rz(-pi) q[2];
rz(-2.0613641) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(-0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85522643) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(-0.65441982) q[0];
x q[1];
rz(-2.084311) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(-1.2531812) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9006173) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(-1.6402628) q[1];
rz(-pi) q[2];
rz(1.0269208) q[3];
sx q[3];
rz(-1.9034991) q[3];
sx q[3];
rz(2.016891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407912) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(-1.1164222) q[2];
sx q[2];
rz(-1.4618256) q[2];
sx q[2];
rz(1.7088919) q[2];
rz(-1.4428044) q[3];
sx q[3];
rz(-0.57590579) q[3];
sx q[3];
rz(-0.87344575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];