OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(-0.85715357) q[0];
sx q[0];
rz(0.13248086) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6462631) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(-2.0908337) q[0];
rz(0.59074596) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(-2.9172446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3092279) q[1];
sx q[1];
rz(-1.698306) q[1];
sx q[1];
rz(-3.068919) q[1];
x q[2];
rz(-1.0971783) q[3];
sx q[3];
rz(-1.3485104) q[3];
sx q[3];
rz(-1.3003295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(1.8923627) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(-0.5805648) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843129) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(2.8312107) q[0];
rz(-pi) q[1];
rz(1.3524019) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(0.44109694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7564056) q[1];
sx q[1];
rz(-2.2573973) q[1];
sx q[1];
rz(2.3910206) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0248236) q[3];
sx q[3];
rz(-1.2516216) q[3];
sx q[3];
rz(1.1899195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(-3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(3.1087648) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2680227) q[0];
sx q[0];
rz(-0.39186726) q[0];
sx q[0];
rz(-0.64143945) q[0];
x q[1];
rz(1.0186586) q[2];
sx q[2];
rz(-1.820192) q[2];
sx q[2];
rz(-2.2160335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9877316) q[1];
sx q[1];
rz(-2.3567833) q[1];
sx q[1];
rz(-1.8042817) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34605108) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(2.3164761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7971928) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(-2.8642505) q[2];
rz(-2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70401496) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(0.2335877) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(2.3707726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17731006) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(1.55127) q[0];
rz(-pi) q[1];
rz(2.0609444) q[2];
sx q[2];
rz(-1.66301) q[2];
sx q[2];
rz(2.8871418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1781436) q[1];
sx q[1];
rz(-1.932486) q[1];
sx q[1];
rz(2.5835035) q[1];
rz(0.869107) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(-1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9757441) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.4404526) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3225587) q[0];
sx q[0];
rz(-1.337262) q[0];
sx q[0];
rz(-2.3755431) q[0];
rz(-2.0427809) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(2.1489378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1483037) q[1];
sx q[1];
rz(-1.9699886) q[1];
sx q[1];
rz(0.94435512) q[1];
x q[2];
rz(0.87336297) q[3];
sx q[3];
rz(-1.0358827) q[3];
sx q[3];
rz(1.8148592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(0.45853841) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(-0.68516723) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72502575) q[0];
sx q[0];
rz(-1.5013114) q[0];
sx q[0];
rz(2.2216703) q[0];
x q[1];
rz(-0.14641996) q[2];
sx q[2];
rz(-0.78260566) q[2];
sx q[2];
rz(2.866982) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.203478) q[1];
sx q[1];
rz(-0.96760975) q[1];
sx q[1];
rz(-0.60738648) q[1];
rz(-pi) q[2];
rz(3.1232883) q[3];
sx q[3];
rz(-2.1566475) q[3];
sx q[3];
rz(2.626112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3926065) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(0.29306456) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.9940631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80619752) q[0];
sx q[0];
rz(-1.1578387) q[0];
sx q[0];
rz(1.6523244) q[0];
x q[1];
rz(0.93716623) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(0.62919754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16312576) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(-0.48467111) q[1];
rz(2.494874) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(-1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-0.914004) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99326837) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-2.3102982) q[0];
rz(1.3759026) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-2.7430699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015394) q[0];
sx q[0];
rz(-1.1050637) q[0];
sx q[0];
rz(-0.17850152) q[0];
x q[1];
rz(-1.6767098) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(0.61818365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9498082) q[1];
sx q[1];
rz(-2.5439918) q[1];
sx q[1];
rz(1.530184) q[1];
x q[2];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(2.2198912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.297696) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(2.0432037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5441355) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(2.068589) q[0];
rz(-0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(-0.84504499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4918412) q[1];
sx q[1];
rz(-1.5713308) q[1];
sx q[1];
rz(0.57926308) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25402756) q[3];
sx q[3];
rz(-1.4048178) q[3];
sx q[3];
rz(2.0167054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(2.6220654) q[2];
rz(1.8064921) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.8241204) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463503) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91613149) q[0];
sx q[0];
rz(-2.7636508) q[0];
sx q[0];
rz(1.4378689) q[0];
x q[1];
rz(-2.9009027) q[2];
sx q[2];
rz(-0.94229892) q[2];
sx q[2];
rz(1.1185874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4911463) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(2.6994929) q[1];
x q[2];
rz(-1.2260776) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(-0.41050875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(1.0591327) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28329904) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(-2.8876866) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(0.72348307) q[2];
sx q[2];
rz(-1.537848) q[2];
sx q[2];
rz(1.6707735) q[2];
rz(-3.0729978) q[3];
sx q[3];
rz(-0.98496901) q[3];
sx q[3];
rz(1.2377644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
