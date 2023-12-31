OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(4.5067956) q[0];
sx q[0];
rz(11.542008) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7919851) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(0.7138568) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65806234) q[2];
sx q[2];
rz(-2.1076638) q[2];
sx q[2];
rz(-1.8368349) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.2092023) q[1];
sx q[1];
rz(-1.1671288) q[1];
rz(-0.23006769) q[3];
sx q[3];
rz(-0.32031968) q[3];
sx q[3];
rz(-0.39428082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(2.1038726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25181928) q[0];
sx q[0];
rz(-2.5368241) q[0];
sx q[0];
rz(-2.7483727) q[0];
rz(-pi) q[1];
rz(2.3677164) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(2.3351923) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2983919) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(-3.1133679) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95275767) q[3];
sx q[3];
rz(-0.85988322) q[3];
sx q[3];
rz(-0.074912138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(2.5684165) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7645435) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(-3.0717875) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8163221) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(-2.6550967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14904505) q[1];
sx q[1];
rz(-0.48679513) q[1];
sx q[1];
rz(2.4942314) q[1];
rz(-2.7335448) q[3];
sx q[3];
rz(-1.9983665) q[3];
sx q[3];
rz(-1.5023295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013214839) q[0];
sx q[0];
rz(-2.2608739) q[0];
sx q[0];
rz(2.8280558) q[0];
rz(-1.5393125) q[2];
sx q[2];
rz(-1.5523947) q[2];
sx q[2];
rz(-3.1061663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4094761) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(-2.9684121) q[1];
rz(-pi) q[2];
rz(-0.94621559) q[3];
sx q[3];
rz(-1.6295625) q[3];
sx q[3];
rz(-1.1138686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5495816) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(0.15079389) q[0];
x q[1];
rz(-2.1826477) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-2.3252955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7534605) q[1];
sx q[1];
rz(-1.7963444) q[1];
sx q[1];
rz(-1.3688341) q[1];
rz(1.0110537) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(-1.884348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(-2.1250336) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(3.0029283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370699) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(-1.5674595) q[0];
rz(-0.20093341) q[2];
sx q[2];
rz(-1.6461419) q[2];
sx q[2];
rz(0.88694015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0181959) q[1];
sx q[1];
rz(-1.0180078) q[1];
sx q[1];
rz(-2.4559896) q[1];
rz(-pi) q[2];
rz(-0.74216446) q[3];
sx q[3];
rz(-1.2483178) q[3];
sx q[3];
rz(3.1165251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(-0.58553186) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666703) q[0];
sx q[0];
rz(-2.2793988) q[0];
sx q[0];
rz(-1.2929582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46632669) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(-1.2517267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1075322) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(2.2570616) q[1];
rz(-2.1950486) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77057225) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(-1.8380941) q[0];
x q[1];
rz(0.89739563) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(0.71133864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1210632) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(1.2709649) q[1];
rz(-pi) q[2];
rz(1.5548607) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006376) q[0];
sx q[0];
rz(-1.3906286) q[0];
sx q[0];
rz(3.0336477) q[0];
rz(1.9753014) q[2];
sx q[2];
rz(-2.1098237) q[2];
sx q[2];
rz(-1.5750969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(0.65111098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3885214) q[3];
sx q[3];
rz(-1.3967447) q[3];
sx q[3];
rz(0.93833246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(-1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-2.5949123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-2.3843117) q[0];
sx q[0];
rz(-1.2454633) q[0];
rz(-pi) q[1];
rz(0.28995138) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(-1.314756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17391275) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(-0.12853865) q[1];
x q[2];
rz(-1.4320847) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(-0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(0.004301087) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.2658723) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(-1.0989582) q[3];
sx q[3];
rz(-1.8660587) q[3];
sx q[3];
rz(2.3793424) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
