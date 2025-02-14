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
rz(-2.9450671) q[0];
sx q[0];
rz(-2.3152469) q[0];
sx q[0];
rz(2.2235121) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(-1.5674051) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188326) q[0];
sx q[0];
rz(-3.1315324) q[0];
sx q[0];
rz(1.2234109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56212933) q[2];
sx q[2];
rz(-1.2527173) q[2];
sx q[2];
rz(-2.0462917) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.955407) q[1];
sx q[1];
rz(-1.8726085) q[1];
sx q[1];
rz(3.0312181) q[1];
rz(1.449578) q[3];
sx q[3];
rz(-1.08076) q[3];
sx q[3];
rz(-3.0090209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7142746) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(1.8753258) q[2];
rz(-1.0424987) q[3];
sx q[3];
rz(-1.5407341) q[3];
sx q[3];
rz(-1.3955759) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365486) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(3.0446766) q[0];
rz(-2.3333683) q[1];
sx q[1];
rz(-2.7327635) q[1];
sx q[1];
rz(1.2867297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42303434) q[0];
sx q[0];
rz(-1.1211532) q[0];
sx q[0];
rz(3.1021523) q[0];
rz(-pi) q[1];
rz(-2.1551844) q[2];
sx q[2];
rz(-0.93743582) q[2];
sx q[2];
rz(1.6478761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5401973) q[1];
sx q[1];
rz(-2.2608065) q[1];
sx q[1];
rz(-2.6371535) q[1];
rz(-1.2169514) q[3];
sx q[3];
rz(-1.8895738) q[3];
sx q[3];
rz(2.6518266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5010995) q[2];
sx q[2];
rz(-0.88397908) q[2];
sx q[2];
rz(-2.7596149) q[2];
rz(-0.52465087) q[3];
sx q[3];
rz(-1.4068406) q[3];
sx q[3];
rz(-2.9978571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3749738) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(2.4554456) q[0];
rz(2.0740017) q[1];
sx q[1];
rz(-2.144404) q[1];
sx q[1];
rz(-0.080726191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2093182) q[0];
sx q[0];
rz(-2.9146412) q[0];
sx q[0];
rz(-1.7448241) q[0];
rz(-pi) q[1];
rz(0.78157921) q[2];
sx q[2];
rz(-0.8197166) q[2];
sx q[2];
rz(2.9224861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9517652) q[1];
sx q[1];
rz(-1.8009342) q[1];
sx q[1];
rz(1.5469502) q[1];
rz(-pi) q[2];
rz(-2.0423546) q[3];
sx q[3];
rz(-2.3009389) q[3];
sx q[3];
rz(1.366758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8778235) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(-0.43886718) q[2];
rz(1.8435439) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(-0.41518655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91442672) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(2.4650204) q[0];
rz(-0.1380955) q[1];
sx q[1];
rz(-2.0510249) q[1];
sx q[1];
rz(2.9409883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56708589) q[0];
sx q[0];
rz(-0.73154035) q[0];
sx q[0];
rz(-1.0163496) q[0];
x q[1];
rz(0.25646957) q[2];
sx q[2];
rz(-0.64531981) q[2];
sx q[2];
rz(0.27094524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.296941) q[1];
sx q[1];
rz(-1.703023) q[1];
sx q[1];
rz(0.6118212) q[1];
rz(-pi) q[2];
rz(-0.225876) q[3];
sx q[3];
rz(-2.568733) q[3];
sx q[3];
rz(1.3536039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(1.7207883) q[2];
rz(3.0270789) q[3];
sx q[3];
rz(-1.4736466) q[3];
sx q[3];
rz(-2.1421471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(0.97212273) q[0];
rz(-0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(1.1824898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33859461) q[0];
sx q[0];
rz(-1.7659682) q[0];
sx q[0];
rz(3.1006569) q[0];
rz(-pi) q[1];
rz(-1.576831) q[2];
sx q[2];
rz(-1.0005992) q[2];
sx q[2];
rz(-1.5046727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69004284) q[1];
sx q[1];
rz(-1.4937972) q[1];
sx q[1];
rz(1.4119215) q[1];
x q[2];
rz(2.413091) q[3];
sx q[3];
rz(-2.430901) q[3];
sx q[3];
rz(-1.817734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9779196) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(-0.53691205) q[2];
rz(2.4162857) q[3];
sx q[3];
rz(-3.091843) q[3];
sx q[3];
rz(-2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86730114) q[0];
sx q[0];
rz(-1.0649571) q[0];
sx q[0];
rz(1.6987479) q[0];
rz(2.9310215) q[1];
sx q[1];
rz(-1.4434394) q[1];
sx q[1];
rz(0.64705667) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4763586) q[0];
sx q[0];
rz(-1.826735) q[0];
sx q[0];
rz(-0.38689918) q[0];
rz(-pi) q[1];
rz(0.55653127) q[2];
sx q[2];
rz(-1.4780296) q[2];
sx q[2];
rz(-1.9829197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3270204) q[1];
sx q[1];
rz(-1.2271735) q[1];
sx q[1];
rz(0.44492857) q[1];
x q[2];
rz(-0.67630597) q[3];
sx q[3];
rz(-2.7190373) q[3];
sx q[3];
rz(0.12862118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0183384) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(1.1798165) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(-3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061851) q[0];
sx q[0];
rz(-1.6949061) q[0];
sx q[0];
rz(0.89757288) q[0];
rz(1.8073742) q[1];
sx q[1];
rz(-1.7709657) q[1];
sx q[1];
rz(-0.99475494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.008381) q[0];
sx q[0];
rz(-1.8803839) q[0];
sx q[0];
rz(0.82010834) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.686196) q[2];
sx q[2];
rz(-2.1242622) q[2];
sx q[2];
rz(1.139251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89284507) q[1];
sx q[1];
rz(-1.1214646) q[1];
sx q[1];
rz(-2.8161418) q[1];
rz(-pi) q[2];
rz(-1.654387) q[3];
sx q[3];
rz(-2.8599742) q[3];
sx q[3];
rz(2.6476988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0576386) q[2];
sx q[2];
rz(-1.3292162) q[2];
sx q[2];
rz(0.29339054) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-0.86881995) q[3];
sx q[3];
rz(-1.4330385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5720125) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(2.8386175) q[0];
rz(-1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(2.4724919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8527864) q[0];
sx q[0];
rz(-1.1256721) q[0];
sx q[0];
rz(1.2034125) q[0];
x q[1];
rz(2.1630493) q[2];
sx q[2];
rz(-2.7251232) q[2];
sx q[2];
rz(-1.5451252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41324612) q[1];
sx q[1];
rz(-1.9921682) q[1];
sx q[1];
rz(2.1284118) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.331393) q[3];
sx q[3];
rz(-1.1576011) q[3];
sx q[3];
rz(1.5511144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9992493) q[2];
sx q[2];
rz(-1.8084904) q[2];
sx q[2];
rz(-0.86714253) q[2];
rz(1.123318) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.860723) q[0];
sx q[0];
rz(-0.46361247) q[0];
sx q[0];
rz(0.11775693) q[0];
rz(-2.2640758) q[1];
sx q[1];
rz(-1.0209457) q[1];
sx q[1];
rz(0.95474517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4539826) q[0];
sx q[0];
rz(-1.1907401) q[0];
sx q[0];
rz(0.65681547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0756726) q[2];
sx q[2];
rz(-1.1469764) q[2];
sx q[2];
rz(1.1413121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2184671) q[1];
sx q[1];
rz(-0.72281853) q[1];
sx q[1];
rz(-1.2413526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8189965) q[3];
sx q[3];
rz(-1.0291417) q[3];
sx q[3];
rz(-1.1198695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2473601) q[2];
sx q[2];
rz(-2.8017513) q[2];
sx q[2];
rz(1.4840688) q[2];
rz(-1.5225211) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(1.9671666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5574944) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(0.75310055) q[0];
rz(1.9901216) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(-1.6023191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26254567) q[0];
sx q[0];
rz(-2.6521054) q[0];
sx q[0];
rz(1.9681843) q[0];
rz(-0.50984211) q[2];
sx q[2];
rz(-0.93889369) q[2];
sx q[2];
rz(-0.84856489) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0286765) q[1];
sx q[1];
rz(-0.43555037) q[1];
sx q[1];
rz(-2.6715379) q[1];
x q[2];
rz(-0.33892858) q[3];
sx q[3];
rz(-2.458161) q[3];
sx q[3];
rz(-3.0290857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3573542) q[2];
sx q[2];
rz(-0.18111649) q[2];
sx q[2];
rz(0.46585807) q[2];
rz(-2.0458131) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(0.54212681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453542) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(-1.0678328) q[1];
sx q[1];
rz(-2.7012431) q[1];
sx q[1];
rz(2.3213097) q[1];
rz(-0.26366641) q[2];
sx q[2];
rz(-1.34524) q[2];
sx q[2];
rz(-3.0312579) q[2];
rz(-1.0505843) q[3];
sx q[3];
rz(-1.0903842) q[3];
sx q[3];
rz(0.36538418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
