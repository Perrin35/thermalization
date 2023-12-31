OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2250741) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(-2.201447) q[0];
rz(-pi) q[1];
rz(0.93809442) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(-3.0245568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7251016) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(-1.7822669) q[1];
rz(-pi) q[2];
rz(1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(-2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(0.73195362) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7322313) q[0];
sx q[0];
rz(-1.1890113) q[0];
sx q[0];
rz(2.2344927) q[0];
rz(-2.5231045) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(2.6478812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6138184) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(2.7223177) q[1];
rz(-1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(1.1211959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7911581) q[0];
sx q[0];
rz(-0.8236304) q[0];
sx q[0];
rz(2.6333548) q[0];
rz(-pi) q[1];
rz(-1.4175225) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(3.0561471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.4157297) q[1];
sx q[1];
rz(-1.9007069) q[1];
sx q[1];
rz(-1.2582474) q[1];
x q[2];
rz(0.49384533) q[3];
sx q[3];
rz(-1.6742799) q[3];
sx q[3];
rz(0.34085694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(-0.63347647) q[0];
x q[1];
rz(2.1999947) q[2];
sx q[2];
rz(-0.59855748) q[2];
sx q[2];
rz(1.4231921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0013106) q[1];
sx q[1];
rz(-1.1661068) q[1];
sx q[1];
rz(2.6005122) q[1];
rz(0.4898407) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7589353) q[0];
sx q[0];
rz(-1.6446582) q[0];
sx q[0];
rz(2.1942684) q[0];
rz(2.2448036) q[2];
sx q[2];
rz(-1.9325581) q[2];
sx q[2];
rz(-2.0828431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7733113) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.7721304) q[1];
x q[2];
rz(0.43907459) q[3];
sx q[3];
rz(-0.77507654) q[3];
sx q[3];
rz(-0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500568) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(0.12775001) q[0];
x q[1];
rz(-0.88271898) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(0.27872745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57950912) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(0.097192055) q[1];
rz(2.0993125) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(-3.022775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(2.8708141) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(2.0628827) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.2111838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-pi) q[1];
rz(1.4326101) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(-1.9463584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83963517) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(-1.1595999) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32254036) q[3];
sx q[3];
rz(-0.99223677) q[3];
sx q[3];
rz(-1.4417779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.221009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85212612) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(2.7702987) q[0];
rz(0.17692716) q[2];
sx q[2];
rz(-1.4104341) q[2];
sx q[2];
rz(-0.71080506) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.358566) q[1];
sx q[1];
rz(-1.6811803) q[1];
sx q[1];
rz(0.13659887) q[1];
x q[2];
rz(-0.82151316) q[3];
sx q[3];
rz(-1.8213846) q[3];
sx q[3];
rz(-0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(3.0158214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98709479) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(2.3214066) q[0];
x q[1];
rz(1.7681098) q[2];
sx q[2];
rz(-1.243405) q[2];
sx q[2];
rz(1.9257853) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51024918) q[1];
sx q[1];
rz(-1.8537087) q[1];
sx q[1];
rz(-0.81088539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5500507) q[3];
sx q[3];
rz(-2.7898443) q[3];
sx q[3];
rz(-1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-2.4323145) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126497) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(1.2055231) q[0];
rz(2.8814949) q[2];
sx q[2];
rz(-2.0195228) q[2];
sx q[2];
rz(-1.5425494) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2503567) q[1];
sx q[1];
rz(-1.6592126) q[1];
sx q[1];
rz(0.18914117) q[1];
x q[2];
rz(-2.8285698) q[3];
sx q[3];
rz(-2.3232197) q[3];
sx q[3];
rz(2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(2.5614212) q[3];
sx q[3];
rz(-2.4709354) q[3];
sx q[3];
rz(1.2448268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
