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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(-1.2641579) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(2.0533419) q[1];
sx q[1];
rz(6.2702141) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91468477) q[0];
sx q[0];
rz(-2.6400262) q[0];
sx q[0];
rz(-0.33059327) q[0];
rz(-pi) q[1];
rz(1.4616724) q[2];
sx q[2];
rz(-1.4089917) q[2];
sx q[2];
rz(1.8617804) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56264191) q[1];
sx q[1];
rz(-1.8899675) q[1];
sx q[1];
rz(2.7719872) q[1];
rz(-pi) q[2];
x q[2];
rz(1.784934) q[3];
sx q[3];
rz(-0.99839568) q[3];
sx q[3];
rz(-0.97808394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.4538572) q[2];
sx q[2];
rz(0.095414735) q[2];
rz(-0.48424193) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(-1.4580844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937623) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(0.65407148) q[0];
rz(-1.3520799) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(-3.0316614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1760565) q[0];
sx q[0];
rz(-1.5586462) q[0];
sx q[0];
rz(-1.5909082) q[0];
rz(-pi) q[1];
rz(-2.1651697) q[2];
sx q[2];
rz(-2.4196673) q[2];
sx q[2];
rz(2.7253369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.667404) q[1];
sx q[1];
rz(-2.0222353) q[1];
sx q[1];
rz(-1.7628481) q[1];
rz(-pi) q[2];
rz(-0.39654332) q[3];
sx q[3];
rz(-1.3351044) q[3];
sx q[3];
rz(-1.9901708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8249417) q[2];
sx q[2];
rz(-1.9330838) q[2];
sx q[2];
rz(1.4671154) q[2];
rz(1.0351099) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255945) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(0.79063928) q[0];
rz(-2.3979893) q[1];
sx q[1];
rz(-1.598282) q[1];
sx q[1];
rz(-0.57949439) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4369219) q[0];
sx q[0];
rz(-1.0365067) q[0];
sx q[0];
rz(0.8922669) q[0];
rz(-pi) q[1];
rz(2.2966301) q[2];
sx q[2];
rz(-0.38533346) q[2];
sx q[2];
rz(0.31491885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.070755006) q[1];
sx q[1];
rz(-2.8808314) q[1];
sx q[1];
rz(1.7160796) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0813339) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(0.34927327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6774595) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(-0.38132384) q[2];
rz(-2.2903806) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(0.73494953) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320936) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(-2.7771948) q[0];
rz(0.20026194) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(-3.0416378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44463377) q[0];
sx q[0];
rz(-1.4958188) q[0];
sx q[0];
rz(3.1211389) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1420629) q[2];
sx q[2];
rz(-1.6661281) q[2];
sx q[2];
rz(0.39458654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32936072) q[1];
sx q[1];
rz(-1.6633004) q[1];
sx q[1];
rz(1.5932371) q[1];
x q[2];
rz(2.3175008) q[3];
sx q[3];
rz(-0.76585117) q[3];
sx q[3];
rz(-1.6834109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.072307) q[2];
sx q[2];
rz(-1.0682718) q[2];
sx q[2];
rz(-3.0266673) q[2];
rz(1.140444) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(0.97525245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9370148) q[0];
sx q[0];
rz(-2.1890722) q[0];
sx q[0];
rz(2.9691147) q[0];
rz(-1.0109488) q[1];
sx q[1];
rz(-2.5806249) q[1];
sx q[1];
rz(0.15484658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5681074) q[0];
sx q[0];
rz(-2.7701039) q[0];
sx q[0];
rz(0.93271352) q[0];
x q[1];
rz(-1.2482951) q[2];
sx q[2];
rz(-2.4716544) q[2];
sx q[2];
rz(-0.33316406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6754969) q[1];
sx q[1];
rz(-1.4267529) q[1];
sx q[1];
rz(2.7269468) q[1];
x q[2];
rz(-2.9633767) q[3];
sx q[3];
rz(-1.6807739) q[3];
sx q[3];
rz(-2.4770155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0871206) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(2.7692914) q[2];
rz(2.7436658) q[3];
sx q[3];
rz(-1.1436661) q[3];
sx q[3];
rz(-0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18846866) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(2.7717223) q[1];
sx q[1];
rz(-1.7431755) q[1];
sx q[1];
rz(-0.63873783) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2681779) q[0];
sx q[0];
rz(-1.8081596) q[0];
sx q[0];
rz(-2.9832212) q[0];
rz(1.9119092) q[2];
sx q[2];
rz(-0.55160597) q[2];
sx q[2];
rz(-0.45845879) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46203278) q[1];
sx q[1];
rz(-2.5278494) q[1];
sx q[1];
rz(-2.2824953) q[1];
rz(-pi) q[2];
rz(-0.79213284) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(1.2611539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0729735) q[2];
sx q[2];
rz(-1.0516473) q[2];
sx q[2];
rz(-0.2529141) q[2];
rz(-0.84826338) q[3];
sx q[3];
rz(-1.6631283) q[3];
sx q[3];
rz(-3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10097583) q[0];
sx q[0];
rz(-2.9315797) q[0];
sx q[0];
rz(0.14895359) q[0];
rz(-1.3111929) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(0.054100903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.220517) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(-1.4709298) q[0];
rz(-0.80246822) q[2];
sx q[2];
rz(-1.5395428) q[2];
sx q[2];
rz(-2.5078338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91840345) q[1];
sx q[1];
rz(-1.4610054) q[1];
sx q[1];
rz(-0.70007433) q[1];
rz(-pi) q[2];
rz(1.6160433) q[3];
sx q[3];
rz(-0.47569573) q[3];
sx q[3];
rz(3.0632116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(-2.1745963) q[2];
rz(-2.4407834) q[3];
sx q[3];
rz(-2.1613224) q[3];
sx q[3];
rz(-0.35082671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0130149) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(-1.5482192) q[0];
rz(-0.37823996) q[1];
sx q[1];
rz(-0.75235569) q[1];
sx q[1];
rz(2.7517448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5704351) q[0];
sx q[0];
rz(-1.3754554) q[0];
sx q[0];
rz(-0.51049149) q[0];
rz(-pi) q[1];
rz(0.016640113) q[2];
sx q[2];
rz(-1.8748054) q[2];
sx q[2];
rz(1.9250411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2521542) q[1];
sx q[1];
rz(-2.5233626) q[1];
sx q[1];
rz(-2.5095224) q[1];
rz(-pi) q[2];
rz(-2.4267219) q[3];
sx q[3];
rz(-0.33239543) q[3];
sx q[3];
rz(-0.22163342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0491911) q[2];
sx q[2];
rz(-2.0202049) q[2];
sx q[2];
rz(1.258705) q[2];
rz(2.8973268) q[3];
sx q[3];
rz(-1.3552856) q[3];
sx q[3];
rz(-2.145483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32023892) q[0];
sx q[0];
rz(-1.9122253) q[0];
sx q[0];
rz(-1.7852596) q[0];
rz(2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(-0.43620268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8863731) q[0];
sx q[0];
rz(-0.56111911) q[0];
sx q[0];
rz(-2.5291689) q[0];
rz(-1.9848083) q[2];
sx q[2];
rz(-1.2591259) q[2];
sx q[2];
rz(3.0729635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.055298565) q[1];
sx q[1];
rz(-2.3196967) q[1];
sx q[1];
rz(2.8248252) q[1];
rz(-pi) q[2];
rz(2.7049191) q[3];
sx q[3];
rz(-1.6323576) q[3];
sx q[3];
rz(1.3832645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.977481) q[2];
sx q[2];
rz(2.5980921) q[2];
rz(0.049526878) q[3];
sx q[3];
rz(-1.4855569) q[3];
sx q[3];
rz(1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.5307584) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(3.0623867) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(2.1150463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3883822) q[0];
sx q[0];
rz(-0.76772707) q[0];
sx q[0];
rz(-1.2305165) q[0];
x q[1];
rz(-1.6658989) q[2];
sx q[2];
rz(-1.9355023) q[2];
sx q[2];
rz(2.022418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0821918) q[1];
sx q[1];
rz(-2.804356) q[1];
sx q[1];
rz(2.853501) q[1];
rz(1.6502569) q[3];
sx q[3];
rz(-1.2209519) q[3];
sx q[3];
rz(-2.3235278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2058699) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(-1.6892461) q[2];
rz(2.3116889) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6185388) q[0];
sx q[0];
rz(-1.9970311) q[0];
sx q[0];
rz(1.507623) q[0];
rz(-1.8745096) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(-2.6020423) q[2];
sx q[2];
rz(-1.0536516) q[2];
sx q[2];
rz(-0.45894844) q[2];
rz(-0.51658705) q[3];
sx q[3];
rz(-2.5174375) q[3];
sx q[3];
rz(0.843912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
