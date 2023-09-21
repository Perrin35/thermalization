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
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(-0.067117604) q[1];
sx q[1];
rz(1.2844515) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0811631) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(-1.1638327) q[0];
rz(-pi) q[1];
rz(-1.60381) q[2];
sx q[2];
rz(-2.0602977) q[2];
sx q[2];
rz(-2.2061493) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40988906) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(1.975592) q[1];
x q[2];
rz(-1.554261) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(-2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54290402) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(2.23404) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.37155) q[2];
sx q[2];
rz(-1.0052048) q[2];
sx q[2];
rz(-2.9233962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-0.60740031) q[1];
sx q[1];
rz(0.022547988) q[1];
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
rz(-pi/2) q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.7012117) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-2.4386141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709135) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(0.12165102) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2146644) q[2];
sx q[2];
rz(-1.0991569) q[2];
sx q[2];
rz(-0.65442649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(-3.0320738) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(-0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(2.3142557) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(2.5783096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-2.1953708) q[0];
sx q[0];
rz(-2.8028691) q[0];
rz(-1.9070542) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(-1.3627571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1174417) q[1];
sx q[1];
rz(-1.5772595) q[1];
sx q[1];
rz(2.8698189) q[1];
rz(1.8654278) q[3];
sx q[3];
rz(-0.8751103) q[3];
sx q[3];
rz(-0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(-2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6974555) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(-0.59209728) q[0];
x q[1];
rz(2.2010872) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(1.3319912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3106766) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(0.78891854) q[1];
x q[2];
rz(-3.1245072) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(2.0320832) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7346142) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(0.91974308) q[0];
rz(-2.6779804) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(-2.6169427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29014153) q[1];
sx q[1];
rz(-2.4212004) q[1];
sx q[1];
rz(3.1300701) q[1];
rz(-pi) q[2];
rz(1.4218876) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(2.7244862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102651) q[0];
sx q[0];
rz(-1.350704) q[0];
sx q[0];
rz(3.1165645) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8646556) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(0.35640946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6228094) q[1];
sx q[1];
rz(-2.2768524) q[1];
sx q[1];
rz(2.519636) q[1];
rz(1.8700637) q[3];
sx q[3];
rz(-0.93942681) q[3];
sx q[3];
rz(-0.33559468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(-3.0292125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880786) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(2.0515576) q[0];
rz(2.8939395) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(2.0695956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1689414) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(2.6886743) q[1];
rz(-pi) q[2];
rz(-1.1376082) q[3];
sx q[3];
rz(-2.8578651) q[3];
sx q[3];
rz(1.197049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(0.91167489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766987) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(2.9093268) q[0];
rz(0.54347221) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(2.6639338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72494353) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(1.3436505) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(-2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(-0.76672673) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.4153597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0534131) q[0];
sx q[0];
rz(-1.3706494) q[0];
sx q[0];
rz(-0.26760177) q[0];
x q[1];
rz(-1.0572817) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(1.2531812) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9006173) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(1.5013298) q[1];
x q[2];
rz(-2.7578027) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407912) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5564556) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.3265058) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(2.1429569) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
